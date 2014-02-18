from pyNastran.f06.f06_formatting import writeFloats13E
from pyNastran.op2.resultObjects.tableObject import RealTableVector, ComplexTableVector, RealTableObject, ComplexTableObject

class RealSPCForcesVector(RealTableVector):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableVector.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n', ]
        #words = ['                                                   V E L O C I T Y   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, pageNum, f)
        return self._write_f06_block(words, header, pageStamp, pageNum, f)


class ComplexSPCForcesVector(ComplexTableVector):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableVector.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                         C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T\n']
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f, is_mag_phase)


class RealSPCForcesObject(RealTableObject):

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, pageNum, f)

        msg = header + words
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            #if not is_all_zeros:
            [dx, dy, dz, rx, ry, rz] = vals2
            msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

        msg.append(pageStamp % pageNum)
        f.write(''.join(msg))
        return pageNum


class ComplexSPCForcesObject(ComplexTableObject):
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                         C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T\n']
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f, is_mag_phase)
