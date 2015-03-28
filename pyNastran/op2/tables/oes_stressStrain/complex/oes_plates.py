from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types
from six.moves import range

from numpy import zeros

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E


class ComplexPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        #self.eType = {}
        #self.element_type = self.data_code['element_type']
        #self.element_name = self.data_code['element_name']

        self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific
        #self.cid = {}  # gridGauss

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self, is_mag_phase):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_nnodes(self):
        if self.element_type in [64, 82, 144]:  # ???, CQUADR, CQUAD4 bilinear
            nnodes = 4 # + 1 centroid
        elif self.element_type in [70, 75]:   #???, CTRIA6
            nnodes = 3 # + 1 centroid
        elif self.element_type in [74, 33]:  # CTRIA3, CQUAD4 linear
            nnodes = 1
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return nnodes + 1  # + 1 centroid

    def build(self):
        #print('data_code = %s' % self.data_code)
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = self.get_nnodes()

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.ntotal = self.nelements * nnodes

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((self.ntotal, 2), 'int32')
        #self.element_cid = zeros((self.nelements, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes * 2 == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)


        self.fiberCurvature = zeros((self.ntotal, 1), 'float32')
        # [oxx, oyy, txy]
        self.data = zeros((self.ntimes, self.ntotal, 3), 'complex64')

    def add_new_eid(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy):
        self.add_eid_sort1(eType, dt, eid, nodeID, fd, oxx, oyy, txy)

    def add(self, dt, eid, gridC, fd, oxx, oyy, txy):
        self.add_eid_sort1(self.element_name, dt, eid, gridC, fd, oxx, oyy, txy)

    def addNewNode(self, dt, eid, gridC, fd, oxx, oyy, txy):
        self.add_eid_sort1(self.element_name, dt, eid, gridC, fd, oxx, oyy, txy)

    def add_eid_sort1(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy):
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s eType=%r dt=%s eid=%s nid=%-5s oxx=%s' % (self.itotal, eType, dt, eid, nodeID, oxx))

        if isinstance(nodeID, string_types):
            nodeID = 0
        self.data[self.itime, self.itotal] = [oxx, oyy, txy]
        self.element_node[self.itotal, :] = [eid, nodeID]  # 0 is center
        self.fiberCurvature[self.itotal] = fd
        #self.ielement += 1
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n' % str(', '.join(self.getHeaders())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.is_stress():
            if self.isFiberDistance():
                gridMsgTemp = ['    ELEMENT              FIBER                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                               '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
                fiberMsgTemp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                                '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
            else:
                gridMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                               '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
                fiberMsgTemp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                                '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
        else:
            if self.isFiberDistance():
                gridMsgTemp = ['    ELEMENT              FIBER                                  - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                               '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
                fiberMsgTemp = ['  ELEMENT       FIBRE                                     - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                                '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
            else:
                gridMsgTemp = ['    ELEMENT              FIBRE                                  - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                               '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
                fiberMsgTemp = ['  ELEMENT       FIBRE                                     - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                                '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']


        if is_mag_phase:
            magReal = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            magReal = ['                                                          (REAL/IMAGINARY)\n', ' \n']

        #if self.isFiberDistance():
            #quadMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -',
                           #'      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY']
        #else:
            #pass

    #'0       100    CEN/8  -2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0'
    #'                       2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0'

        #quadMsg  = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'] + formWord + quadMsgTemp
        #quadrMsg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )'] + formWord + quadMsgTemp
        #quad8Msg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )'] + formWord + quadMsgTemp

        ## TODO: validation on header formatting...
        if self.is_stress():
            cquad4_bilinear = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n']
            cquad4_linear = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
            cquad8 = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
            cquadr = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
            ctria3 = ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']  # good
            ctria6 = ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
            ctriar = ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']
        else:
            cquad4_bilinear = ['                C O M P L E X   S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n']
            cquad4_linear = ['                C O M P L E X   S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']
            cquad8 = ['                C O M P L E X   S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
            cquadr = ['                C O M P L E X   S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
            ctria3 = ['                C O M P L E X   S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']
            ctria6 = ['                C O M P L E X   S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
            ctriar = ['                C O M P L E X   S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']

        msg = []
        is_bilinear = False
        if self.element_name == 'CQUAD4':
            if self.element_type == 144:
                is_bilinear = True
                msg += cquad4_linear + magReal + gridMsgTemp
            elif self.element_type == 33:
                is_bilinear = False
                msg += cquad4_bilinear + magReal + fiberMsgTemp
        elif self.element_name == 'CQUAD8':
            msg += cquad8 + magReal + gridMsgTemp
            is_bilinear = True
        elif self.element_name == 'CQUADR':
            msg += cquadr + magReal + gridMsgTemp
            is_bilinear = True

        elif self.element_name == 'CTRIA3':
            msg += ctria3 + fiberMsgTemp
        elif self.element_name == 'CTRIA6':
            msg += ctria6 + gridMsgTemp
            is_bilinear = True
        elif self.element_name == 'CTRIAR':
            msg += ctriar + gridMsgTemp
            is_bilinear = True
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

        nnodes = self.get_nnodes()
        return msg, nnodes, is_bilinear

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        (msg_temp, nnodes, is_bilinear) = self.get_f06_header(is_mag_phase)

        (ntimes, ntotal, six) = self.data.shape
        for itime in range(ntimes):
            dt = self._times[itime]

            #print('eids=', eids)

            dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dtLine
            msg = header + msg_temp
            f.write('\n'.join(msg))

            if self.element_type == 144: # CQUAD4 bilinear
                self.writeF06_Quad4_BilinearTransient(f, itime, 4, is_mag_phase, 'CEN/4')
            elif self.get_element_type == 33:  # CQUAD4 linear
                self.writeF06_Tri3Transient(f, itime, 4, is_mag_phase, 'CEN/4')
            elif self.element_type == 74: # CTRIA3
                self.writeF06_Tri3Transient(f, itime, 3, is_mag_phase, 'CEN/3')
            elif self.element_type == 64:  #CQUAD8
                self.writeF06_Quad4_BilinearTransient(f, itime, 5, is_mag_phase, 'CEN/8')
            elif self.element_type == 82:  # CQUADR
                self.writeF06_Quad4_BilinearTransient(f, itime, 5, is_mag_phase, 'CEN/8')
            elif self.element_type ==  70:  # CTRIAR
                self.writeF06_Quad4_BilinearTransient(f, itime, 3, is_mag_phase, 'CEN/3')
            elif self.element_type == 75:  # CTRIA6
                self.writeF06_Quad4_BilinearTransient(f, itime, 3, is_mag_phase, 'CEN/6')
            else:
                raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def writeF06_Tri3Transient(self, f, itime, n, is_magnitude_phase, cen):
        """
        CQUAD4 linear
        CTRIA3
        """
        fds = self.fiberCurvature[:, 0]
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        ilayer0 = True
        for eid, node, fd, doxx, doyy, dtxy in zip(eids, nodes, fds, oxx, oyy, txy):
            vals, is_all_zeros = writeFloats13E([fd])
            fd = vals[0]
            ([oxxr, oyyr, txyr,
              oxxi, oyyi, txyi,], is_all_zeros) = writeImagFloats13E([doxx, doyy, dtxy], is_magnitude_phase)

            if ilayer0:    # TODO: assuming 2 layers?
                f.write('0  %6i   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % (eid, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            else:
                f.write('   %6s   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % ('', fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            ilayer0 = not ilayer0

    def writeF06_Quad4_BilinearTransient(self, f, itime, n, is_magnitude_phase, cen):
        """
        CQUAD4 bilinear
        CQUAD8
        CTRIAR
        CTRIA6
        """
        fds = self.fiberCurvature[:, 0]
        oxx = self.data[itime, :, 0]
        oyy = self.data[itime, :, 1]
        txy = self.data[itime, :, 2]

        eids = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        ilayer0 = True
        for eid, node, fd, doxx, doyy, dtxy in zip(eids, nodes, fds, oxx, oyy, txy):
            vals, is_all_zeros = writeFloats13E([fd])
            fd = vals[0]
            ([oxxr, oyyr, txyr,
              oxxi, oyyi, txyi,], is_all_zeros) = writeImagFloats13E([doxx, doyy, dtxy], is_magnitude_phase)

            if node == 0 and ilayer0:
                f.write('0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n' % (eid, cen, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            elif ilayer0:    # TODO: assuming 2 layers?
                f.write('   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n' % ('', node, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            else:
                f.write('   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n\n' % ('', '', fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            ilayer0 = not ilayer0


class ComplexPlateStressArray(ComplexPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    #def is_stress(self):
        #return True

    #def is_strain(self):
        #return False

    def getHeaders(self):
        return ['oxx', 'oyy', 'txy']

class ComplexPlateStrainArray(ComplexPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def getHeaders(self):
        return ['exx', 'eyy', 'exy']

    #def is_stress(self):
        #return False

    #def is_strain(self):
        #return True


class ComplexPlateStress(StressObject):
    """
    ::

                  C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )
                                                            (REAL/IMAGINARY)

      ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -
        ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY
  0       100    CEN/8  -2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0
                         2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}

        #self.append_data_member('sCodes','s_code')
        #print "self.s_code = ",self.s_code
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.fiberCurvature = {}
        self.oxx = {}
        self.oyy = {}
        self.txy = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                self.addNewNode = self.addNewNodeSort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2
            self.addNewNode = self.addNewNodeSort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, fiberCurvature, oxx, oyy, txy\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType == 'CTRIA3':
                    (eType, eid, f1, ox1, oy1, txy1,
                     f2, ox2, oy2, txy2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'CEN/3': [f1, f2]}
                    self.oxx[eid] = {'CEN/3': [ox1, ox2]}
                    self.oyy[eid] = {'CEN/3': [oy1, oy2]}
                    self.txy[eid] = {'CEN/3': [txy1, txy2]}
                elif eType == 'CQUAD4':
                    if len(line) == 19:  # Centroid - bilinear
                        (eType, eid, nid, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.oxx[eid] = {nid: [ox1, ox2]}
                        self.oyy[eid] = {nid: [oy1, oy2]}
                        self.txy[eid] = {nid: [txy1, txy2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        nid = 'CEN/4'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.oxx[eid] = {nid: [ox1, ox2]}
                        self.oyy[eid] = {nid: [oy1, oy2]}
                        self.txy[eid] = {nid: [txy1, txy2]}
                    elif len(line) == 17:  # Bilinear
                        (nid, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        self.fiberCurvature[eid][nid] = [f1, f2]
                        self.oxx[eid][nid] = [ox1, ox2]
                        self.oyy[eid][nid] = [oy1, oy2]
                        self.txy[eid][nid] = [txy1, txy2]
                    else:
                        assert len(line) == 19, 'len(line)=%s' % len(line)
                        raise NotImplementedError()
                else:
                    raise NotImplementedError('line=%s not supported...' % line)
            return
        raise NotImplementedError('transient results not supported')

    def delete_transient(self, dt):
        #del self.fiberCurvature[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.txy[dt]

    def get_transients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        #self.fiberCurvature[dt] = {}
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.txy[dt] = {}

    def add_new_eid_sort1(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy):
        #msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" %(dt,eid,nodeID,fd,oxx,oyy,txy)
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%s" % (
            dt, eid, nodeID, fd, oxx)
        #print msg
        #if eid in self.oxx[dt]:
        #    return self.add(eid,nodeID,fd,oxx,oyy,txy)

        if 0: # this is caused by superelements
            if dt in self.oxx and eid in self.oxx[dt]:  # SOL200, erase the old result
                #nid = nodeID
                #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s" %(dt,eid,nodeID,str(self.fiberCurvature[eid][nid]),str(self.oxx[dt][eid][nid])))
                self.delete_transient(dt)
                self.add_new_transient(dt)

        assert isinstance(eid, int)
        assert isinstance(nodeID, int), nodeID
        self.eType[eid] = eType
        if dt not in self.oxx:
            self.add_new_transient(dt)
        self.fiberCurvature[eid] = {nodeID: [fd]}
        self.oxx[dt][eid] = {nodeID: [oxx]}
        self.oyy[dt][eid] = {nodeID: [oyy]}
        self.txy[dt][eid] = {nodeID: [txy]}
        #if nodeID == 0:
            #raise ValueError(msg)

    def add_sort1(self, dt, eid, nodeID, fd, oxx, oyy, txy):
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" % (
            dt, eid, nodeID, fd, oxx, oyy, txy)
        assert eid is not None
        assert isinstance(nodeID, int), nodeID
        self.fiberCurvature[eid][nodeID].append(fd)
        self.oxx[dt][eid][nodeID].append(oxx)
        self.oyy[dt][eid][nodeID].append(oyy)
        self.txy[dt][eid][nodeID].append(txy)
        #if nodeID == 0:
            #raise ValueError(msg)

    def addNewNodeSort1(self, dt, eid, nodeID, fd, oxx, oyy, txy):
        assert eid is not None
        assert isinstance(nodeID, int), nodeID
        #msg = "eid=%s nodeID=%s fd=%g oxx=%s oyy=%s txy=%s" % (
            #eid, nodeID, fd, oxx, oyy, txy)
        #print msg
        #assert nodeID not in self.oxx[dt][eid]
        self.fiberCurvature[eid][nodeID] = [fd]
        self.oxx[dt][eid][nodeID] = [oxx]
        self.oyy[dt][eid][nodeID] = [oyy]
        self.txy[dt][eid][nodeID] = [txy]
        #if nodeID == 0:
            #raise ValueError(msg)

    def getHeaders(self):
        if self.is_fiber_distance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']
        headers += ['oxx', 'oyy', 'txy']
        return headers

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        assert f is not None
        if len(self.eType) == 0:
            raise RuntimeError('The result object is empty')
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase)

        if is_mag_phase:
            magReal = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            magReal = ['                                                          (REAL/IMAGINARY)\n', ' \n']

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -',
                           '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY']
        else:
            pass

#'0       100    CEN/8  -2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0'
#'                       2.500000E-02    0.0          /  0.0             0.0          /  0.0             0.0          /  0.0'
        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None

        quadMsg  = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'] + formWord + quadMsgTemp
        quadrMsg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )'] + formWord + quadMsgTemp
        quad8Msg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )'] + formWord + quadMsgTemp

        eTypes = self.eType.values()
        msg_packs = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        valid_types = ['CTRIA3', 'CTRIA6', 'CTRIAR',
                       'CQUAD4', 'CQUAD8', 'CQUADR']
        etypes, ordered_etypes = self.getOrderedETypes(valid_types)

        msg = []
        for eType in etypes:
            eids = ordered_etypes[eType]
            if eids:
                eids.sort()
                msg_pack = msg_packs[eType]

                msg += header + msg_pack
                if eType in ['CQUAD4']:
                    if is_bilinear:
                        for eid in eids:
                            out = self.writeF06_Quad4_Bilinear(eid, 4, is_mag_phase, 'CEN/4')
                            msg.append(out)
                    else:
                        for eid in eids:
                            out = self.writeF06_Tri3(eid, 4, is_mag_phase)
                            msg.append(out)
                elif eType in ['CTRIA3']:
                    for eid in eids:
                        out = self.writeF06_Tri3(eid, 3, is_mag_phase)
                        msg.append(out)
                elif eType in ['CQUAD8']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 8, is_mag_phase, 'CEN/8')
                        msg.append(out)
                elif eType in ['CTRIAR']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 3, is_mag_phase, 'CEN/3')
                        msg.append(out)
                elif eType in ['CTRIA6']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 3, is_mag_phase, 'CEN/6')
                        msg.append(out)
                else:
                    raise NotImplementedError('eType = %r' % eType)  # CQUAD8, CTRIA6

                msg.append(page_stamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        return page_num - 1

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        if is_mag_phase:
            magReal = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            magReal = ['                                                          (REAL/IMAGINARY)\n \n']

        if self.isFiberDistance():
            gridMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                           '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiberMsgTemp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                            '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
        else:
            gridMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                           '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiberMsgTemp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                            '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']

        #if self.isFiberDistance():
            #quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
            #               '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            #triMsgTemp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
            #              '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
        #else:
            #quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
            #               '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            #triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
            #              '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]
        #quadMsg  = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )']+magReal+fiberMsgTemp
        #quadrMsg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )']+magReal+fiberMsgTemp
        #quad8Msg = ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )']+magReal+fiberMsgTemp
        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None

        eTypes = self.eType.values()
        dts = self.oxx.keys()
        dt = dts[0]
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[dt][kkey].keys()
            is_bilinear = True
            quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + magReal + gridMsgTemp
            if len(ekey) == 1:
                is_bilinear = False
                quadMsg = header + ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + magReal + fiberMsgTemp

        if 'CQUAD8' in eTypes:
            quad8Msg = header + ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + magReal + gridMsgTemp

        if 'CQUADR' in eTypes:
            quadrMsg = header + ['                C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + magReal + gridMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + magReal + fiberMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + magReal + gridMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + magReal + gridMsgTemp

        msg_packs = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (eTypes, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        dts = self.oxx.keys()
        dts.sort()
        for eType in eTypes:
            eids = orderedETypes[eType]
            if eids:
                msg_pack = msg_packs[eType]
                eids.sort()
                if eType in ['CQUAD4']:
                    if is_bilinear:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.data_code[
                                'name'], dt)
                            msg += header + msg_pack
                            for eid in eids:
                                out = self.writeF06_Quad4_BilinearTransient(dt, eid, 4, is_mag_phase, 'CEN/4')
                                msg.append(out)
                            msg.append(page_stamp % page_num)
                            page_num += 1
                            f.write(''.join(msg))
                            msg = ['']
                    else:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                            msg += header + msg_pack
                            for eid in eids:
                                out = self.writeF06_Tri3Transient(dt, eid, 4, is_mag_phase)
                                msg.append(out)
                            msg.append(page_stamp % page_num)
                            page_num += 1
                            f.write(''.join(msg))
                            msg = ['']
                elif eType in ['CTRIA3']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                        msg += header + msg_pack
                        for eid in eids:
                            out = self.writeF06_Tri3Transient(dt, eid, 3, is_mag_phase)
                            msg.append(out)
                        msg.append(page_stamp % page_num)
                        page_num += 1
                        f.write(''.join(msg))
                        msg = ['']
                elif eType in ['CTRIAR']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        msg += header + msg_pack
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 3, is_mag_phase, 'CEN/3')
                            msg.append(out)
                        msg.append(page_stamp % page_num)
                        page_num += 1
                        f.write(''.join(msg))
                        msg = ['']
                elif eType in ['CTRIA6']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        msg += header + msg_pack
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 3, is_mag_phase, 'CEN/6')
                            msg.append(out)
                        msg.append(page_stamp % page_num)
                        page_num += 1
                        f.write(''.join(msg))
                        msg = ['']
                elif eType in ['CQUAD8']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        msg += header + msg_pack
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 8, is_mag_phase, 'CEN/8')
                            msg.append(out)
                        msg.append(page_stamp % page_num)
                        page_num += 1
                        f.write(''.join(msg))
                        msg = ['']
                elif eType in ['CQUADR']:
                    if is_bilinear:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.data_code[
                                'name'], dt)
                            msg += header + msg_pack
                            for eid in eids:
                                out = self.writeF06_Quad4_BilinearTransient(dt, eid, 4, is_mag_phase, 'CEN/4')
                                msg.append(out)
                            msg.append(page_stamp % page_num)
                            page_num += 1
                            f.write(''.join(msg))
                            msg = ['']
                    else:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                            msg += header + msg_pack
                            for eid in eids:
                                out = self.writeF06_Tri3Transient(dt, eid, 4, is_mag_phase)
                                msg.append(out)
                            msg.append(page_stamp % page_num)
                            page_num += 1
                            f.write(''.join(msg))
                            msg = ['']
                else:
                    # CQUAD8, CTRIA6
                    raise NotImplementedError('eType = %r' % eType)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        return page_num - 1

    def writeF06_Quad4_Bilinear(self, eid, n, is_mag_phase, cen):
        msg = ''
        #k = self.oxx[eid].keys()
        #k.remove(cen)
        #k.sort()
        #nids = [cen] + k
        nids = sorted(self.oxx[eid].keys())
        for nid in nids:
            for ilayer in range(len(self.oxx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][ilayer]
                oxx = self.oxx[eid][nid][ilayer]
                oyy = self.oyy[eid][nid][ilayer]
                txy = self.txy[eid][nid][ilayer]
                ([fd, oxxr, oyyr, txyr,
                  fdi, oxxi, oyyi, txyi], is_all_zeros) = writeImagFloats13E([fd, oxx, oyy, txy], is_mag_phase)

                if nid == cen and ilayer == 0:
                    msg += '0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % (eid, cen, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                elif ilayer == 0:
                    msg += '   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % ('', nid, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                elif ilayer == 1:
                    msg += '   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n\n' % ('', '', fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                else:
                    raise Exception('Invalid option for cquad4')
        return msg

    def writeF06_Quad4_BilinearTransient(self, dt, eid, n, is_mag_phase, cen):
        msg = ''
        #k = self.oxx[dt][eid].keys()
        ##cen = 'CEN/' + str(n)
        #k.remove(cen)
        #k.sort()
        #nids = [cen] + k
        nids = sorted(self.oxx[dt][eid].keys())
        for nid in nids:
            for ilayer in range(len(self.oxx[dt][eid][nid])):
                fdr = self.fiberCurvature[eid][nid][ilayer]
                oxx = self.oxx[dt][eid][nid][ilayer]
                oyy = self.oyy[dt][eid][nid][ilayer]
                txy = self.txy[dt][eid][nid][ilayer]
                ([fdr, oxxr, oyyr, txyr,
                  fdi, oxxi, oyyi, txyi], is_all_zeros) = writeImagFloats13E([fd, oxx, oyy, txy], is_mag_phase)

                if nid == cen and ilayer == 0:
                    msg += '0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % (eid, cen, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                elif ilayer == 0:
                    msg += '   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % ('', nid, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                elif ilayer == 1:
                    msg += '   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n\n' % ('', '', fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                else:
                    #msg += '   %8s %8s  %13E  %13E %13E %13E   %8.4F  %13E %13E %13E\n' %('','',  fd,oxx,oyy,txy)
                    raise Exception('Invalid option for cquad4')
        return msg

    def writeF06_Tri3(self, eid, n, is_mag_phase):
        msg = ''
        #k = self.oxx[eid].keys()
        #cen = 'CEN/' + str(n)
        #k.remove(cen)
        #k.sort()
        #nids = [cen] + k
        nids = sorted(self.oxx[eid].keys())
        for nid in nids:
            for ilayer in range(len(self.oxx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][ilayer]
                oxx = self.oxx[eid][nid][ilayer]
                oyy = self.oyy[eid][nid][ilayer]
                txy = self.txy[eid][nid][ilayer]
                ([fd, oxx, oyy, txy], is_all_zeros) = writeFloats13E([fd, oxx, oyy, txy])

                if ilayer == 0:
                    # TODO: how is this valid?
                    msg += '0G  %6i   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %s\n' % (eid, fd, oxx, oyy, txy)
                else:
                    msg += ' H  %6s   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %s\n' % ('', fd, oxx, oyy, txy)
        return msg

    def writeF06_Tri3Transient(self, dt, eid, n, is_mag_phase):
        msg = ''
        #k = self.oxx[dt][eid].keys()
        #cen = 'CEN/' + str(n)
        #k.remove(cen)
        #k.sort()
        #nids = [cen] + k
        nids = sorted(self.oxx[dt][eid].keys())
        for nid in nids:
            for ilayer in range(len(self.oxx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][ilayer]
                oxx = self.oxx[dt][eid][nid][ilayer]
                oyy = self.oyy[dt][eid][nid][ilayer]
                txy = self.txy[dt][eid][nid][ilayer]
                ([fd, oxxr, oyyr, txyr,
                  fdi, oxxi, oyyi, txyi], is_all_zeros) = writeImagFloats13E([fd, oxx, oyy, txy], is_mag_phase)

                if ilayer == 0:
                    msg += '0  %6i   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % (eid, fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                else:
                    msg += '   %6s   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % ('', fd, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
        return msg


class ComplexPlateStrain(StrainObject):
    """
    ::

      # ??? - is this just 11
      ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)
        ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES

      # s_code=11
                             S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)
        ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       VON MISES

      # s_code=15
                             S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
      ELEMENT      FIBER                STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)
        ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES

      # s_code=10
                             S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)          MAX
        ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.fiberCurvature = {}
        self.exx = {}
        self.eyy = {}
        self.exy = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                self.addNewNode = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2
            self.addNewNode = self.addNewNodeSort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.exx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, fiberCurvature, exx, eyy, exy\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType == 'CTRIA3':
                    (eType, eid, f1, ex1, ey1, exy1,
                     f2, ex2, ey2, exy2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'CEN/3': [f1, f2]}
                    self.exx[eid] = {'CEN/3': [ex1, ex2]}
                    self.eyy[eid] = {'CEN/3': [ey1, ey2]}
                    self.exy[eid] = {'CEN/3': [exy1, exy2]}
                elif eType == 'CQUAD4':
                    if len(line) == 19:  # Centroid - bilinear
                        (eType, eid, nid, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.exx[eid] = {nid: [ex1, ex2]}
                        self.eyy[eid] = {nid: [ey1, ey2]}
                        self.exy[eid] = {nid: [exy1, exy2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        nid = 'CEN/4'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.exx[eid] = {nid: [ex1, ex2]}
                        self.eyy[eid] = {nid: [ey1, ey2]}
                        self.exy[eid] = {nid: [exy1, exy2]}
                    elif len(line) == 17:  # Bilinear node
                        (nid, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        self.fiberCurvature[eid][nid] = [f1, f2]
                        self.exx[eid][nid] = [ex1, ex2]
                        self.eyy[eid][nid] = [ey1, ey2]
                        self.txy[eid][nid] = [exy1, exy2]
                    else:
                        assert len(line) == 19, 'len(line)=%s' % len(line)
                        raise NotImplementedError()
                else:
                    raise NotImplementedError('line=%s not supported...' % line)
            return
        raise NotImplementedError('transient results not supported')

    def delete_transient(self, dt):
        #del self.fiberCurvature[dt]
        del self.exx[dt]
        del self.eyy[dt]
        del self.exy[dt]

    def get_transients(self):
        k = self.exx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #self.fiberCurvature = {}
        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.exy[dt] = {}

    def add_new_eid_sort1(self, eType, dt, eid, nodeID, curvature, exx, eyy, exy):
        msg = "eid=%s nodeID=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, nodeID, curvature, exx, eyy, exy)

        #if nodeID != 'C':  # centroid
            #assert 0 < nodeID < 1000000000, 'nodeID=%s %s' % (nodeID, msg)

        if 1: # this is caused by superelements
            if dt in self.exx and eid in self.exx[dt]:  # SOL200, erase the old result
                #nid = nodeID
                #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s" %(dt,eid,nodeID,str(self.fiberCurvature[eid][nid]),str(self.exx[dt][eid][nid])))
                self.delete_transient(dt)
                self.add_new_transient(dt)

        #if eid in self.exx[dt]:  # SOL200, erase the old result
            #nid = nodeID
            #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s" %(dt,eid,nodeID,str(self.fiberCurvature[eid][nid]),str(self.oxx[dt][eid][nid]))
            #self.delete_transient(dt)
            #self.add_new_transient()

        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[dt][eid] = {nodeID: [exx]}
        self.eyy[dt][eid] = {nodeID: [eyy]}
        self.exy[dt][eid] = {nodeID: [exy]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add_sort1(self, dt, eid, nodeID, curvature, exx, eyy, exy):
        msg = "eid=%s nodeID=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, nodeID, curvature, exx, eyy, exy)
        #if nodeID != 'C':  # centroid
            #assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)

        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[dt][eid][nodeID].append(exx)
        self.eyy[dt][eid][nodeID].append(eyy)
        self.exy[dt][eid][nodeID].append(exy)
        if nodeID == 0:
            raise ValueError(msg)


    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']

        headers += ['exx', 'eyy', 'exy']
        if self.isVonMises():
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        assert f is not None
        if len(self.eType) == 0:
            raise RuntimeError('The result object is empty')
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f)
        raise RuntimeError('this can never happen')

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER                STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]

        quadMsg = None
        quad8Msg = None
        quadrMsg = None
        triMsg = None
        tri6Msg = None
        trirMsg = None

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            ElemKey = eTypes.index('CQUAD4')
            #print qkey
            eid = self.eType.keys()[ElemKey]
            #print "self.oxx = ",self.oxx
            #print "eid=%s" %(eid)
            dt = self.exx.keys()[0]
            #print "dt=%s" %(dt)
            nLayers = len(self.exx[dt][eid])
            #print "elementKeys = ",elementKeys
            is_bilinear = True
            quadMsg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if nLayers == 1:
                is_bilinear = False
                quadMsg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        msg = []
        msg_packs = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        dts = self.exx.keys()
        dts.sort()
        for eType in typesOut:
            eids = orderedETypes[eType]
            if eids:
                eids.sort()
                eType = self.eType[eid]
                if eType in ['CQUAD4']:
                    if is_bilinear:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                            for eid in eids:
                                out = self.writeF06_Quad4_BilinearTransient(dt, eid, 4, is_mag_phase, 'CEN/4')
                                msg.append(out)
                            msg.append(page_stamp % page_num)
                            f.write(''.join(msg))
                            msg = ['']
                            page_num += 1
                    else:
                        for dt in dts:
                            header[1] = ' %s = %10.4E\n' % (self.data_code[
                                'name'], dt)
                            for eid in eids:
                                out = self.writeF06_Tri3Transient(dt, eid, 4, is_mag_phase)
                                msg.append(out)
                            msg.append(page_stamp % page_num)
                            f.write(''.join(msg))
                            msg = ['']
                            page_num += 1
                elif eType in ['CTRIA3']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        for eid in eids:
                            out = self.writeF06_Tri3Transient(dt, eid, 3, is_mag_phase)
                            msg.append(out)
                        msg.append(page_stamp % page_num)
                        f.write(''.join(msg))
                        msg = ['']
                        page_num += 1
                elif eType in ['CQUAD8']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 5, is_mag_phase, 'CEN/8')
                            msg.append(out)
                        msg.append(page_stamp % page_num)
                        f.write(''.join(msg))
                        msg = ['']
                        page_num += 1
                elif eType in ['CTRIAR']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 3, is_mag_phase, 'CEN/3')
                            msg.append(out)
                        msg.append(page_stamp % page_num)
                        f.write(''.join(msg))
                        msg = ['']
                        page_num += 1
                elif eType in ['CTRIA6']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (
                            self.data_code['name'], dt)
                        for eid in eids:
                            out = self.writeF06_Quad4_BilinearTransient(dt, eid, 3, is_mag_phase, 'CEN/6')
                            msg.append(out)
                        msg.append(page_stamp % page_num)
                        f.write(''.join(msg))
                        msg = ['']
                        page_num += 1
                else:
                    raise NotImplementedError('eType = %r' % eType)  # CQUAD8, CTRIA6
        return page_num - 1

    def writeF06_Quad4_BilinearTransient(self, dt, eid, n, is_mag_phase, cen):
        msg = ''
        #k = self.exx[dt][eid].keys()
        #k.remove(cen)
        #k.sort()
        #nids = [cen] + k
        nids = sorted(self.exx[dt][eid].keys())
        for nid in nids:
            for ilayer in range(len(self.exx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][ilayer]
                exx = self.exx[dt][eid][nid][ilayer]
                eyy = self.eyy[dt][eid][nid][ilayer]
                exy = self.exy[dt][eid][nid][ilayer]
                ([fd, exxr, eyyr, exyr,
                  fdi, exxi, eyyi, exyi], is_all_zeros) = writeImagFloats13E([fd, exx, eyy, exy], is_mag_phase)

                if nid == cen and ilayer == 0:
                    msg += '0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % (eid, cen, fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
                elif ilayer == 0:
                    msg += '   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % ('', nid, fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
                elif ilayer == 1:
                    msg += '   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n\n' % ('', '', fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
                else:
                    raise Exception('Invalid option for cquad4')
        return msg

    def writeF06_Tri3Transient(self, dt, eid, n, is_mag_phase):
        msg = ''
        #k = self.exx[dt][eid].keys()
        #cen = 'CEN/' + str(n)
        #k.remove(cen)
        #k.sort()
        #nids = [cen] + k
        nids = sorted(self.exx[dt][eid].keys())
        for nid in nids:
            for ilayer in range(len(self.exx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][ilayer]
                exx = self.exx[dt][eid][nid][ilayer]
                eyy = self.eyy[dt][eid][nid][ilayer]
                exy = self.exy[dt][eid][nid][ilayer]

                ([fd, exxr, eyyr, exyr,
                  fdi, exxi, eyyi, exyi], is_all_zeros) = writeImagFloats13E([fd, exx, eyy, exy], is_mag_phase)
                if ilayer == 0:
                    msg += '0  %6i   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % (eid, fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
                else:
                    msg += '   %6s   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % ('', fd, exxr, exxi, eyyr, eyyi, exyr, exyi)
        return msg

