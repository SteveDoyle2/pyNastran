from numpy import array

from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E


class EigenVectorObject(TableObject):  # approach_code=2, sort_code=0, thermal=0
    """
    ::

      EIGENVALUE =  6.158494E+07
          CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
             1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
          2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
          2003      G      0.0            0.0            0.0            0.0            0.0            0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, imode, read_mode):
        TableObject.__init__(self, data_code, is_sort1, isubcase, imode, read_mode)
        self.gridTypes = {}

    def read_f06_data(self, data_code, data):
        imode = data_code['mode']
        if imode not in self.translations:
            self.update_mode(data_code, imode)

        for line in data:
            (nid, gridType, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[nid] = gridType
            self.translations[imode][nid] = array([t1, t2, t3])
            self.rotations[imode][nid] = array([r1, r2, r3])
        assert self.eigrs[-1] == data_code['eigr']

    def eigenvalues(self):
        return self.eigrs

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'eigenvectors'
        return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        """
        ::

          EIGENVALUE =  6.158494E+07
              CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

          POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
                 1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
              2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
              2003      G      0.0            0.0            0.0            0.0            0.0            0.0
        """
        msg = []
        hasCycle = hasattr(self, 'mode_cycle')

        #print "self.eigrs =", self.eigrs
        #print "dir",dir(self)
        nmodes, nnodes, modes = self._get_shape()
        #ntotal = nmodes * nnodes
        data = self.data
        #mode0 = modes[0]
        #data.ix[mode0]
        ndata = len(self.data)
        i = 0
        while i < ndata:
            index = self.data.index[i]
            (mode, node_id) = index
            imode = modes.index(mode)
            msg += header
            freq = self.eigrs[imode]
            msg.append('%16s = %13E\n' % ('EIGENVALUE', freq))

            if hasCycle:
                msg.append('%16s = %13E         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, imode))
            else:
                msg.append('                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % imode)

            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')

            gridType = 'G'
            #gridType = self.gridTypes[nodeID]
            #ix = data.index[]
            mode_old = mode
            while mode == mode_old:
                index = self.data.index[i]
                (mode, node_id) = index
                #nodeID = data[index]['node_id']
                dx = data['T1'][index]
                dy = data['T2'][index]
                dz = data['T3'][index]
                rx = data['R1'][index]
                ry = data['R2'][index]
                rz = data['R3'][index]

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (node_id, gridType, dx, dy, dz, rx, ry, rz.rstrip()))
                i += 1
                try:
                    index = self.data.index[i]
                except:
                    break
                (mode, node_id) = index
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1


class RealEigenVectorObject(ScalarObject):  # approach_code=2, sort_code=0, thermal=0
    """
    ::

                                           R E A L   E I G E N V E C T O R   N O .          1
        POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
               1      G      0.0            0.0            0.0            0.0            1.260264E-01   0.0
               7      G      9.999849E-01   0.0            6.728968E-03   0.0            8.021386E-03   0.0
    """
    def __init__(self, data_code, isubcase, imode, read_mode):
        self.shape = {}
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.caseVal = mode
        #print "mode = %s" % imode
        self.caseVal = self.getUnsteadyValue()
        self.set_data_members()
        self.gridTypes = {}

    def modes(self):
        return sorted(self.translations.keys())

    def eigenvalues(self):
        return self.eigrs

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'eigenvectors'
        return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        """
        ::

          EIGENVALUE =  6.158494E+07
                                             R E A L   E I G E N V E C T O R   N O .          1

          POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
                 1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
              2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
              2003      G      0.0            0.0            0.0            0.0            0.0            0.0
        """
        msg = []
        #print self.data_code
        for i, (iMode, eigVals) in enumerate(sorted(self.translations.iteritems())):
            msg += header
            freq = self.eigrs[i]
            msg.append('%16s = %12E\n' % ('EIGENVALUE', freq))
            msg.append('                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' % (iMode))
            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            for nodeID, translation in sorted(eigVals.iteritems()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dx, dy, dz, rx, ry, rz.rstrip()))
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        return (''.join(msg), pageNum - 1)


class ComplexEigenVectorObject(ComplexTableObject):  # approach_code=2, sort_code=0, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, imode, read_mode):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, imode, read_mode)
        self.shape = {}

    def eigenvalues(self):
        return sorted(self.translations.keys())

    def write_matlab(self, isubcase, f=None, is_mag_phase=False):
        name = 'complex_eigenvectors'
        return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        msg = []
        #print self.data_code
        hasCycle = hasattr(self, 'mode_cycle')
        for i, (iMode, eigVals) in enumerate(sorted(self.translations.iteritems())):
            msg += header
            freq = self.eigrs[i]
            #freq = 0.0
            msg.append('%16s = %12E\n' % ('EIGENVALUE', freq))
            if hasCycle:
                msg.append('%16s = %12E          C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, iMode))
            else:
                msg.append('                                         C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % (iMode))
            msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
            for nodeID, displacement in sorted(eigVals.iteritems()):
                rotation = self.rotations[iMode][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx, dy, dz) = displacement
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = writeImagFloats13E(vals, is_mag_phase)
                [dxr, dyr, dzr, rxr, ryr, rzr, dxi, dyi,
                    dzi, rxi, ryi, rzi] = vals2

                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, gridType, dxr, dyr, dzr, rxr, ryr, rzr.rstrip()))
                msg.append('%14s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % ('', '', dxi, dyi, dzi, rxi, ryi, rzi.rstrip()))

            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1