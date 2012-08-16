import sys
#import copy

# pyNastran
from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject


class AccelerationObject(TableObject):  # approachCode=11, thermal=0
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        TableObject.__init__(self, dataCode, isSort1, iSubcase, dt)

    def writeMatlab(self, iSubcase, f=None, isMagPhase=False):
        name = 'accelerations'
        if self.nonlinearFactor is None:
            return self._writeMatlab(name, iSubcase, f)
        else:
            return self._writeMatlabTransient(name, iSubcase, f)

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.getTableMarker()
        return self._writeF06Block(words, header, pageStamp, pageNum, f)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.getTableMarker()
        return self._writeF06TransientBlock(words, header, pageStamp, pageNum, f)

    def __repr__(self):
        #return ''
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---ACCELERATIONS---\n'
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
        msg = '---TRANSIENT ACCELERATIONS---\n'
        msg += self.writeHeader()

        for dt, translations in sorted(self.translations.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
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
                    ###
                msg += '\n'
            ###
        return msg


class ComplexAccelerationObject(ComplexTableObject):  # tableCode=11, approachCode=???
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        ComplexTableObject.__init__(self, dataCode, isSort1, iSubcase, dt)

    def writeMatlab(self, iSubcase, f=None, isMagPhase=False):
        name = 'accelerations'
        if self.nonlinearFactor is None:
            return self._writeMatlab(name, iSubcase, f, isMagPhase)
        else:
            return self._writeMatlabTransient(name, iSubcase, f, isMagPhase)

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f, isMagPhase)

        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._writeF06Block(words, header, pageStamp, pageNum, f, isMagPhase)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._writeF06TransientBlock(words, header, pageStamp, pageNum, f, isMagPhase)

    def __repr__(self):
        return self.writeF06(['', '', ''], 'PAGE ', 1)[0]
        #if self.nonlinearFactor is not None:
            #return self.__reprTransient__()

        msg = '---COMPLEX ACCELERATIONS---\n'
        #if self.nonlinearFactor is not None:
        #    msg += '%s = %g\n' %(self.dataCode['name'],self.dt)
        headers = ['DxReal', 'DxImag', 'DyReal', 'DyImag', 'DzReal', 'DyImag', 'RxReal', 'RxImag', 'RyReal', 'RyImag', 'RzReal', 'RzImag']
        msg += '%-10s ' % ('nodeID')
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for freq, translations in sorted(self.translations.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], freq)

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
