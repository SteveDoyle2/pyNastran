from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types
from six.moves import range

from numpy import zeros

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E, get_key0


class ComplexPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

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

    def get_nnodes(self):
        return get_nnodes(self)

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

    def add_new_eid_sort1(self, eType, dt, eid, node_id, fdr, oxx, oyy, txy):
        self.add_eid_sort1(eType, dt, eid, node_id, fdr, oxx, oyy, txy)

    def add_sort1(self, dt, eid, gridC, fdr, oxx, oyy, txy):
        self.add_eid_sort1(self.element_name, dt, eid, gridC, fdr, oxx, oyy, txy)

    def add_new_node_sort1(self, dt, eid, gridc, fdr, oxx, oyy, txy):
        self.add_eid_sort1(self.element_name, dt, eid, gridc, fdr, oxx, oyy, txy)

    def add_eid_sort1(self, eType, dt, eid, node_id, fdr, oxx, oyy, txy):
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #print('itotal=%s eType=%r dt=%s eid=%s nid=%-5s oxx=%s' % (self.itotal, eType, dt, eid, node_id, oxx))

        assert isinstance(node_id, int), node_id
        self.data[self.itime, self.itotal] = [oxx, oyy, txy]
        self.element_node[self.itotal, :] = [eid, node_id]  # 0 is center
        self.fiberCurvature[self.itotal] = fdr
        #self.ielement += 1
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
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
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg_temp, nnodes, is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1)

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]

            dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dtLine
            msg = header + msg_temp
            f.write('\n'.join(msg))

            if self.element_type == 144: # CQUAD4 bilinear
                self._write_f06_quad4_bilinear_transient(f, itime, 4, is_mag_phase, 'CEN/4')
            elif self.element_type == 33:  # CQUAD4 linear
                self._write_f06_tri3_transient(f, itime, 4, is_mag_phase, 'CEN/4')
            elif self.element_type == 74: # CTRIA3
                self._write_f06_tri3_transient(f, itime, 3, is_mag_phase, 'CEN/3')
            elif self.element_type == 64:  #CQUAD8
                self._write_f06_quad4_bilinear_transient(f, itime, 5, is_mag_phase, 'CEN/8')
            elif self.element_type == 82:  # CQUADR
                self._write_f06_quad4_bilinear_transient(f, itime, 5, is_mag_phase, 'CEN/8')
            elif self.element_type == 70:  # CTRIAR
                self._write_f06_quad4_bilinear_transient(f, itime, 3, is_mag_phase, 'CEN/3')
            elif self.element_type == 75:  # CTRIA6
                self._write_f06_quad4_bilinear_transient(f, itime, 3, is_mag_phase, 'CEN/6')
            else:
                raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def _write_f06_tri3_transient(self, f, itime, n, is_magnitude_phase, cen):
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
        for eid, node, fdr, doxx, doyy, dtxy in zip(eids, nodes, fds, oxx, oyy, txy):
            vals, is_all_zeros = writeFloats13E([fdr])
            fdr = vals[0]
            ([oxxr, oyyr, txyr,
              oxxi, oyyi, txyi,], is_all_zeros) = writeImagFloats13E([doxx, doyy, dtxy], is_magnitude_phase)

            if ilayer0:    # TODO: assuming 2 layers?
                f.write('0  %6i   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % (eid, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            else:
                f.write('   %6s   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % ('', fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            ilayer0 = not ilayer0

    def _write_f06_quad4_bilinear_transient(self, f, itime, n, is_magnitude_phase, cen):
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
            fdr = vals[0]
            ([oxxr, oyyr, txyr,
              oxxi, oyyi, txyi,], is_all_zeros) = writeImagFloats13E([doxx, doyy, dtxy], is_magnitude_phase)

            if node == 0 and ilayer0:
                f.write('0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n' % (eid, cen, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            elif ilayer0:    # TODO: assuming 2 layers?
                f.write('   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n' % ('', node, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            else:
                f.write('   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s / %s\n\n' % ('', '', fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi))
            ilayer0 = not ilayer0

def _get_plate_msg(self, is_mag_phase=True, is_sort1=True):
    if self.is_von_mises():
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    if self.is_stress():
        if self.is_fiber_distance():
            grid_msg_temp = ['    ELEMENT              FIBER                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                             '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiber_msg_temp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                              '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
        else:
            grid_msg_temp = ['    ELEMENT              FIBRE                                  - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                             '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiber_msg_temp = ['  ELEMENT       FIBRE                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -\n',
                              '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
    else:
        if self.is_fiber_distance():
            grid_msg_temp = ['    ELEMENT              FIBER                                  - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                             '      ID      GRID-ID   DISTANCE                 NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiber_msg_temp = ['  ELEMENT       FIBRE                                     - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                              '    ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']
        else:
            grid_msg_temp = ['    ELEMENT              FIBRE                                  - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                             '      ID      GRID-ID   CURVATURE                NORMAL-X                        NORMAL-Y                       SHEAR-XY\n']
            fiber_msg_temp = ['  ELEMENT       FIBRE                                     - STRAINS IN ELEMENT  COORDINATE SYSTEM -\n',
                              '    ID.       CURVATURE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY\n']


    if is_mag_phase:
        mag_real = ['                                                         (MAGNITUDE/PHASE)\n \n']
    else:
        mag_real = ['                                                          (REAL/IMAGINARY)\n', ' \n']

    #if self.is_fiber_distance():
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
    if self.element_type == 144: # CQUAD4
        is_bilinear = True
        msg += cquad4_linear + mag_real + grid_msg_temp
    elif self.element_type == 33: # CQUAD4
        is_bilinear = False
        msg += cquad4_bilinear + mag_real + fiber_msg_temp
    elif self.element_type == 64:  #CQUAD8
        msg += cquad8 + mag_real + grid_msg_temp
        is_bilinear = True
    elif self.element_type == 82:  # CQUADR
        msg += cquadr + mag_real + grid_msg_temp
        is_bilinear = True

    elif self.element_type == 74: # CTRIA3
        msg += ctria3 + fiber_msg_temp
    elif self.element_type == 75:  # CTRIA6
        msg += ctria6 + grid_msg_temp
        is_bilinear = True
    elif self.element_type == 70:  # CTRIAR
        msg += ctriar + grid_msg_temp
        is_bilinear = True
    else:
        raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

    nnodes = get_nnodes(self)
    return msg, nnodes, is_bilinear

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

class ComplexPlateStressArray(ComplexPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        return ['oxx', 'oyy', 'txy']

class ComplexPlateStrainArray(ComplexPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        return ['exx', 'eyy', 'exy']


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
                    cen = 0 # CEN/3
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {cen : [f1, f2]}
                    self.oxx[eid] = {cen : [ox1, ox2]}
                    self.oyy[eid] = {cen : [oy1, oy2]}
                    self.txy[eid] = {cen : [txy1, txy2]}
                elif eType == 'CQUAD4':
                    if len(line) == 19:  # Centroid - bilinear
                        (eType, eid, node_id, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        self.eType[eid] = eType
                        assert isinstance(node_id, int), node_id
                        self.fiberCurvature[eid] = {node_id : [f1, f2]}
                        self.oxx[eid] = {node_id : [ox1, ox2]}
                        self.oyy[eid] = {node_id : [oy1, oy2]}
                        self.txy[eid] = {node_id : [txy1, txy2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        node_id = 0 # CEN/4
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {node_id : [f1, f2]}
                        self.oxx[eid] = {node_id : [ox1, ox2]}
                        self.oyy[eid] = {node_id : [oy1, oy2]}
                        self.txy[eid] = {node_id : [txy1, txy2]}
                    elif len(line) == 17:  # Bilinear
                        (node_id, f1, ox1, oy1, txy1,
                         f2, ox2, oy2, txy2) = line
                        assert isinstance(node_id, int), node_id
                        self.fiberCurvature[eid][node_id] = [f1, f2]
                        self.oxx[eid][node_id] = [ox1, ox2]
                        self.oyy[eid][node_id] = [oy1, oy2]
                        self.txy[eid][node_id] = [txy1, txy2]
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

    def add_new_eid_sort1(self, eType, dt, eid, node_id, fdr, oxx, oyy, txy):
        #msg = "dt=%s eid=%s node_id=%s fdr=%g oxx=%s oyy=%s txy=%s" % (dt, eid, node_id, fdr, oxx, oyy, txy)
        #msg = "dt=%s eid=%s node_id=%s fdr=%g oxx=%s" % (dt, eid, node_id, fdr, oxx)
        #if eid in self.oxx[dt]:
        #    return self.add(eid, node_id, fdr, oxx, oyy, txy)

        assert isinstance(eid, int)
        assert isinstance(node_id, int), node_id
        self.eType[eid] = eType
        if dt not in self.oxx:
            self.add_new_transient(dt)
        self.fiberCurvature[eid] = {node_id : [fdr]}
        self.oxx[dt][eid] = {node_id : [oxx]}
        self.oyy[dt][eid] = {node_id : [oyy]}
        self.txy[dt][eid] = {node_id : [txy]}

    def add_new_node_sort1(self, dt, eid, node_id, fdr, oxx, oyy, txy):
        msg = "eid=%s node_id=%s fdr=%g oxx=%s oyy=%s txy=%s" % (eid, node_id, fdr, oxx, oyy, txy)
        assert eid is not None, eid
        assert isinstance(node_id, int), node_id
        assert node_id > -1, msg
        #assert node_id not in self.oxx[dt][eid]
        self.fiberCurvature[eid][node_id] = [fdr]
        self.oxx[dt][eid][node_id] = [oxx]
        self.oyy[dt][eid][node_id] = [oyy]
        self.txy[dt][eid][node_id] = [txy]

    def add_sort1(self, dt, eid, node_id, fdr, oxx, oyy, txy):
        msg = "dt=%s eid=%s node_id=%s fdr=%g oxx=%s oyy=%s txy=%s" % (dt, eid, node_id, fdr, oxx, oyy, txy)
        assert eid is not None, eid
        assert isinstance(node_id, int), node_id
        assert node_id > -1, msg

        self.fiberCurvature[eid][node_id].append(fdr)
        self.oxx[dt][eid][node_id].append(oxx)
        self.oyy[dt][eid][node_id].append(oyy)
        self.txy[dt][eid][node_id].append(txy)

    def _get_headers(self):
        if self.is_fiber_distance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']
        headers += ['oxx', 'oyy', 'txy']
        return headers

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        assert f is not None
        if len(self.eType) == 0:
            raise RuntimeError('The result object is empty')
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        raise RuntimeError('this can never happen')

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg_pack, nnodes, is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1)

        msg = []
        dts = list(self.oxx.keys())
        dts.sort()
        dt0 = dts[0]
        dt_name = self.data_code['name']

        eids = sorted(self.oxx[dt0].keys())
        eids.sort()

        if self.element_type == 144: # CQUAD4
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (dt_name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_quad4_bilinear_transient(dt, eid, 4, is_mag_phase, 'CEN/4')
                    msg.append(out)
                msg.append(page_stamp % page_num)
                page_num += 1
                f.write(''.join(msg))
                msg = ['']
        elif self.element_type == 33: # CQUAD4
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (dt_name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_tri3_transient(dt, eid, 4, is_mag_phase)
                    msg.append(out)
                msg.append(page_stamp % page_num)
                page_num += 1
                f.write(''.join(msg))
                msg = ['']
        elif self.element_type == 74: # CTRIA3
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (dt_name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_tri3_transient(dt, eid, 3, is_mag_phase)
                    msg.append(out)
                msg.append(page_stamp % page_num)
                page_num += 1
                f.write(''.join(msg))
                msg = ['']
        elif self.element_type == 70:  # CTRIAR
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (dt_name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_quad4_bilinear_transient(dt, eid, 3, is_mag_phase, 'CEN/3')
                    msg.append(out)
                msg.append(page_stamp % page_num)
                page_num += 1
                f.write(''.join(msg))
                msg = ['']
        elif self.element_type == 75:  # CTRIA6
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (dt_name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_quad4_bilinear_transient(dt, eid, 3, is_mag_phase, 'CEN/6')
                    msg.append(out)
                msg.append(page_stamp % page_num)
                page_num += 1
                f.write(''.join(msg))
                msg = ['']
        elif self.element_type == 64:  #CQUAD8
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (dt_name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_quad4_bilinear_transient(dt, eid, 8, is_mag_phase, 'CEN/8')
                    msg.append(out)
                msg.append(page_stamp % page_num)
                page_num += 1
                f.write(''.join(msg))
                msg = ['']
        elif self.element_type == 82:  # CQUADR
            if is_bilinear:
                for dt in dts:
                    header[1] = ' %s = %10.4E\n' % (dt_name, dt)
                    msg += header + msg_pack
                    for eid in eids:
                        out = self._write_f06_quad4_bilinear_transient(dt, eid, 4, is_mag_phase, 'CEN/4')
                        msg.append(out)
                    msg.append(page_stamp % page_num)
                    page_num += 1
                    f.write(''.join(msg))
                    msg = ['']
            else:
                for dt in dts:
                    header[1] = ' %s = %10.4E\n' % (dt_name, dt)
                    msg += header + msg_pack
                    for eid in eids:
                        out = self._write_f06_tri3_transient(dt, eid, 4, is_mag_phase)
                        msg.append(out)
                    msg.append(page_stamp % page_num)
                    page_num += 1
                    f.write(''.join(msg))
                    msg = ['']
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return page_num - 1

    def _write_f06_quad4_bilinear_transient(self, dt, eid, n, is_mag_phase, cen):
        """
        CQUAD4 bilinear
        CQUAD8
        CTRIAR
        CTRIA6
        """
        msg = ''
        nids = sorted(self.oxx[dt][eid].keys())
        for node_id in nids:
            for ilayer in range(len(self.oxx[dt][eid][node_id])):
                fdr = self.fiberCurvature[eid][node_id][ilayer]
                oxx = self.oxx[dt][eid][node_id][ilayer]
                oyy = self.oyy[dt][eid][node_id][ilayer]
                txy = self.txy[dt][eid][node_id][ilayer]
                ([fdr, oxxr, oyyr, txyr,
                  fdi, oxxi, oyyi, txyi], is_all_zeros) = writeImagFloats13E([fdr, oxx, oyy, txy], is_mag_phase)

                if node_id == 0 and ilayer == 0:
                    msg += '0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % (eid, cen, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                elif ilayer == 0:
                    msg += '   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % ('', node_id, fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                elif ilayer == 1:
                    msg += '   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n\n' % ('', '', fdr, oxxr, oxxi, oyyr, oyyi, txyr, txyi)
                else:
                    #msg += '   %8s %8s  %13E  %13E %13E %13E   %8.4F  %13E %13E %13E\n' % ('', '',  fd, oxx, oyy, txy)
                    raise RuntimeError('Invalid option for cquad4')
        return msg

    def _write_f06_tri3_transient(self, dt, eid, n, is_mag_phase):
        msg = ''
        nids = sorted(self.oxx[dt][eid].keys())
        for node_id in nids:
            for ilayer in range(len(self.oxx[dt][eid][node_id])):
                fd = self.fiberCurvature[eid][node_id][ilayer]
                oxx = self.oxx[dt][eid][node_id][ilayer]
                oyy = self.oyy[dt][eid][node_id][ilayer]
                txy = self.txy[dt][eid][node_id][ilayer]
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
                    cen = 0 # CEN/3
                    self.fiberCurvature[eid] = {cen : [f1, f2]}
                    self.exx[eid] = {cen : [ex1, ex2]}
                    self.eyy[eid] = {cen : [ey1, ey2]}
                    self.exy[eid] = {cen : [exy1, exy2]}
                elif eType == 'CQUAD4':
                    if len(line) == 19:  # Centroid - bilinear
                        (eType, eid, node_id, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        assert isinstance(node_id, int), node_id
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {node_id : [f1, f2]}
                        self.exx[eid] = {node_id : [ex1, ex2]}
                        self.eyy[eid] = {node_id : [ey1, ey2]}
                        self.exy[eid] = {node_id : [exy1, exy2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        node_id = 0 # CEN/4
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {node_id : [f1, f2]}
                        self.exx[eid] = {node_id : [ex1, ex2]}
                        self.eyy[eid] = {node_id : [ey1, ey2]}
                        self.exy[eid] = {node_id : [exy1, exy2]}
                    elif len(line) == 17:  # Bilinear node
                        (node_id, f1, ex1, ey1, exy1,
                         f2, ex2, ey2, exy2) = line
                        assert isinstance(node_id, int), node_id
                        self.fiberCurvature[eid][node_id] = [f1, f2]
                        self.exx[eid][node_id] = [ex1, ex2]
                        self.eyy[eid][node_id] = [ey1, ey2]
                        self.exy[eid][node_id] = [exy1, exy2]
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

    def add_new_eid_sort1(self, eType, dt, eid, node_id, curvature, exx, eyy, exy):
        msg = "eid=%s node_id=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, node_id, curvature, exx, eyy, exy)

        #if node_id != 'C':  # centroid
            #assert 0 < node_id < 1000000000, 'node_id=%s %s' % (node_id, msg)
        assert isinstance(node_id, int), node_id
        assert node_id > -1, msg

        if dt not in self.exx:
            msg = "dt=%s eid=%s node_id=%s" % (dt, eid, node_id)
            self.add_new_transient(dt)

        #if eid in self.exx[dt]:  # SOL200, erase the old result
            #nid = node_id
            #msg = "dt=%s eid=%s node_id=%s fd=%s oxx=%s" %(dt,eid,node_id,str(self.fiberCurvature[eid][node_id]),str(self.oxx[dt][eid][node_id]))
            #self.delete_transient(dt)
            #self.add_new_transient()

        self.eType[eid] = eType
        self.add_new_node_sort1(dt, eid, node_id, curvature, exx, eyy, exy)

    def add_new_node_sort1(self, dt, eid, node_id, curvature, exx, eyy, exy):
        self.fiberCurvature[eid] = {node_id : [curvature]}
        self.exx[dt][eid] = {node_id : [exx]}
        self.eyy[dt][eid] = {node_id : [eyy]}
        self.exy[dt][eid] = {node_id : [exy]}

    def add_sort1(self, dt, eid, node_id, curvature, exx, eyy, exy):
        msg = "eid=%s node_id=%s curvature=%g exx=%s eyy=%s exy=%s" % (
            eid, node_id, curvature, exx, eyy, exy)
        assert isinstance(node_id, int), node_id
        assert node_id > -1, msg

        self.fiberCurvature[eid][node_id].append(curvature)
        self.exx[dt][eid][node_id].append(exx)
        self.eyy[dt][eid][node_id].append(eyy)
        self.exy[dt][eid][node_id].append(exy)


    def _get_headers(self):
        if self.is_fiber_distance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']

        headers += ['exx', 'eyy', 'exy']
        if self.is_von_mises():
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        assert f is not None
        if len(self.eType) == 0:
            raise RuntimeError('The result object is empty')
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        raise RuntimeError('this can never happen')

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg_pack, nnodes, is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1)

        msg = []
        dts = list(self.exx.keys())
        dts.sort()
        ntimes = len(dts)
        dt0 = dts[0]
        eids = sorted(self.exx[dt0])

        name = self.data_code['name']
        if self.element_type == 144: # CQUAD4 bilinear
            if is_bilinear:
                for dt in dts:
                    header[1] = ' %s = %10.4E\n' % (name, dt)
                    msg += header + msg_pack
                    for eid in eids:
                        out = self._write_f06_quad4_bilinear_transient(dt, eid, 4, is_mag_phase, 'CEN/4')
                        msg.append(out)
                    msg.append(page_stamp % page_num)
                    f.write(''.join(msg))
                    msg = ['']
                    page_num += 1
        elif self.element_type == 33:  # CQUAD4 linear
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_tri3_transient(dt, eid, 4, is_mag_phase)
                    msg.append(out)
                msg.append(page_stamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        elif self.element_type == 74: # CTRIA3
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_tri3_transient(dt, eid, 3, is_mag_phase)
                    msg.append(out)
                msg.append(page_stamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        elif self.element_type == 64:  #CQUAD8
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_quad4_bilinear_transient(dt, eid, 5, is_mag_phase, 'CEN/8')
                    msg.append(out)
                msg.append(page_stamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        elif self.element_type == 70:  # CTRIAR
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_quad4_bilinear_transient(dt, eid, 3, is_mag_phase, 'CEN/3')
                    msg.append(out)
                msg.append(page_stamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        elif self.element_type == 75:  # CTRIA6
            for dt in dts:
                header[1] = ' %s = %10.4E\n' % (name, dt)
                msg += header + msg_pack
                for eid in eids:
                    out = self._write_f06_quad4_bilinear_transient(dt, eid, 3, is_mag_phase, 'CEN/6')
                    msg.append(out)
                msg.append(page_stamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        else:
    #FREQUENCY =  1.800000E+01
        #C O M P L E X     S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
            #(MAGNITUDE/PHASE)

    #ELEMENT       FIBER                                     -  STRAINS IN ELEMENT  COORDINATE SYSTEM -
        #ID.        DISTANCE                  NORMAL-X                          NORMAL-Y                         SHEAR-XY
#0   51001   -2.380000E-03      2.471426E-08 /   4.2486           2.662035E-08 / 183.7563           5.579315E-07 /   3.6197
    #2.380000E-03      1.048796E-07 /   3.7203           1.379607E-08 / 183.9262           5.700372E-07 /   3.6174

            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return page_num - 1

    def _write_f06_quad4_bilinear_transient(self, dt, eid, n, is_mag_phase, cen):
        """
        CQUAD4 bilinear
        CQUAD8
        CTRIAR
        CTRIA6
        """
        msg = ''
        nids = sorted(self.exx[dt][eid].keys())
        for node_id in nids:
            for ilayer in range(len(self.exx[dt][eid][node_id])):
                fdr = self.fiberCurvature[eid][node_id][ilayer]
                exx = self.exx[dt][eid][node_id][ilayer]
                eyy = self.eyy[dt][eid][node_id][ilayer]
                exy = self.exy[dt][eid][node_id][ilayer]
                ([fdr, exxr, eyyr, exyr,
                  fdi, exxi, eyyi, exyi], is_all_zeros) = writeImagFloats13E([fdr, exx, eyy, exy], is_mag_phase)

                if node_id == 0 and ilayer == 0:
                    msg += '0  %8i %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % (eid, cen, fdr, exxr, exxi, eyyr, eyyi, exyr, exyi)
                elif ilayer == 0:
                    msg += '   %8s %8i  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n' % ('', node_id, fdr, exxr, exxi, eyyr, eyyi, exyr, exyi)
                elif ilayer == 1:
                    msg += '   %8s %8s  %-13s   %-13s / %-13s   %-13s / %-13s   %-13s /   %s\n\n' % ('', '', fdr, exxr, exxi, eyyr, eyyi, exyr, exyi)
                else:
                    raise RuntimeError('Invalid option for cquad4')
        return msg

    def _write_f06_tri3_transient(self, dt, eid, n, is_mag_phase):
        msg = ''
        nids = sorted(self.exx[dt][eid].keys())
        for node_id in nids:
            for ilayer in range(len(self.exx[dt][eid][node_id])):
                fdr = self.fiberCurvature[eid][node_id][ilayer]
                exx = self.exx[dt][eid][node_id][ilayer]
                eyy = self.eyy[dt][eid][node_id][ilayer]
                exy = self.exy[dt][eid][node_id][ilayer]

                ([fdr, exxr, eyyr, exyr,
                  fdi, exxi, eyyi, exyi], is_all_zeros) = writeImagFloats13E([fdr, exx, eyy, exy], is_mag_phase)
                if ilayer == 0:
                    msg += '0  %6i   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % (eid, fdr, exxr, exxi, eyyr, eyyi, exyr, exyi)
                else:
                    msg += '   %6s   %-13s     %-13s / %-13s     %-13s / %-13s     %-13s / %s\n' % ('', fdr, exxr, exxi, eyyr, eyyi, exyr, exyi)
        return msg
