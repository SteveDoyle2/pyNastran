#pylint disable=C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types
from six.moves import zip, range
from itertools import count
from numpy import zeros, searchsorted, ravel

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, writeFloats8p4F

class RealPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                self.addNewNode = self.addNewNodeSort1
        else:
            raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2
            #self.addNewNode = self.addNewNodeSort2

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def is_bilinear(self):
        if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
            return False
        elif self.element_type in [144]:  # CQUAD4
            return False
        raise NotImplementedError(self.element_type)

    def build(self):
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 74:
            nnodes_per_element = 1
        elif self.element_type == 144:
            nnodes_per_element = 5

        self.nnodes = nnodes_per_element
        #self.nelements //= nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.ntotal, 2), dtype='int32')

        #[fd, oxx, oyy, txy, angle, majorP, minorP, ovm]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def add_new_eid(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        self.add_new_eid_sort1(eType, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)

    def add_new_eid_sort1(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print("------------------")
        #print("isStress =", self.isStress())
        #msg = "i=%s eType=%s dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (
            #self.itotal,  eType, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print(msg)

        assert isinstance(eid, int)
        self._times[self.itime] = dt
        assert isinstance(nodeID, string_types), nodeID
        if isinstance(nodeID, string_types):
            nodeID = 0
        self.element_node[self.itotal, :] = [eid, nodeID]
        self.data[self.itime, self.itotal, :] = [fd, oxx, oyy, txy, angle, majorP, minorP, ovm]
        self.itotal += 1
        self.ielement += 1

    def addNewNode(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        if isinstance(nodeID, string_types):
            nodeID = 0
        self.add_sort1(dt, eid, nodeID, fd, oxx, oyy, txy, angle,
                            majorP, minorP, ovm)


    def add(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        if isinstance(nodeID, string_types):
            nodeID = 0
        self.add_sort1(dt, eid, nodeID, fd, oxx, oyy, txy, angle,
                            majorP, minorP, ovm)

    def addNewNodeSort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        if isinstance(nodeID, string_types):
            nodeID = 0
        self.add_sort1(dt, eid, nodeID, fd, oxx, oyy, txy, angle,
                            majorP, minorP, ovm)

    def add_sort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print('addNewNodeSort1')
        assert eid is not None
        msg = "i=%s dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (
            self.itotal, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print(msg)
        if isinstance(nodeID, string_types):
            nodeID = 0
        #assert isinstance(nodeID, int), nodeID
        self.element_node[self.itotal, :] = [eid, nodeID]
        self.data[self.itime, self.itotal, :] = [fd, oxx, oyy, txy, angle, majorP, minorP, ovm]
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.nnodes
        ntotal = self.ntotal
        nlayers = 2
        nelements = self.ntotal // self.nnodes // 2

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, nlayers, ntotal))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n'
                       % (self.__class__.__name__, nelements, nnodes, nlayers, ntotal))
            ntimes_word = 1
        #msg.append('  data.shape=%s' % str(self.data.shape))
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element types: %s\n  ' % ', '.join(self.element_names))
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        ctria3_msg, ctria6_msg, cquad4_msg, cquad8_msg = self._get_msgs()
        if 'CTRIA3' in self.element_name and self.element_type == 74:
            msg = ctria3_msg
            nnodes = 3
        elif 'CQUAD4' in self.element_name and self.element_type == 33:
            msg = cquad4_msg
            nnodes = 4
        elif 'CTRIA6' in self.element_name and self.element_type == 0:
            msg = ctria6_msg
            nnodes = 6
        elif 'CQUAD8' in self.element_name and self.element_type == 0:
            msg = cquad8_msg
            nnodes = 8
            raise RuntimeError('can these be bilinear???')
        else:
            raise NotImplementedError(self.element_name)

        return self.element_name, nnodes, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        msg, nnodes, is_bilinear = self._get_msgs()

        # write the f06
        (ntimes, ntotal, four) = self.data.shape

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        cen_word = 'CEN/%i' % nnodes
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            if self.nonlinear_factor is not None:
                dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                header[1] = dtLine
                if hasattr(self, 'eigr'):
                    header[2] = ' %14s = %12.6E\n' % ('EIGENVALUE', self.eigrs[itime])
            f.write(''.join(header + msg))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[fd, oxx, oyy, txy, angle, majorP, minorP, ovm]
            fd = self.data[itime, :, 0]
            oxx = self.data[itime, :, 1]
            oyy = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            angle = self.data[itime, :, 4]
            majorP = self.data[itime, :, 5]
            minorP = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]

            # loop over all the elements
            for i, eid, nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi in zip(
                count(), eids, nids, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
                ([fdi, oxxi, oyyi, txyi, major, minor, ovmi],
                 is_all_zeros) = writeFloats13E([fdi, oxxi, oyyi, txyi, major, minor, ovmi])
                #f.write([eidi, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi])
                iLayer = i % 2
                # tria3
                if self.element_type == 74:
                    if iLayer == 0:
                        f.write('0  %6i   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (eid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    else:
                        f.write('   %6s   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % ('', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))

                elif self.element_type == 144:
                    # bilinear
                    if nid == 0 and iLayer == 0:  # CEN
                        f.write('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (eid, cen_word, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif iLayer == 0:
                        f.write('   %8s %8i  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % ('', nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif iLayer == 1:
                        f.write('   %8s %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n\n' % ('', '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                else:
                    raise RuntimeError(self.element_type)

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealPlateStressArray(RealPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def isStress(self):
        return True
    def isStrain(self):
        return False

    def get_headers(self):
        if self.isFiberDistance():
            fd = 'fiber_distance'
        else:
            fd = 'fiber_curvature'

        if self.isVonMises():
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = [fd, 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', ovm]
        return headers

    def _get_msgs(self):
        if self.isVonMises():
            von_mises = 'VON MISES'
        else:
            von_mises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        else:
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                           '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]

        #=============================

        isBilinear = False
        cquad4_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp
        cquad8_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        ## TODO: can cquad8s be bilinear???
        isBilinear = True
        cquadr_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
        cquad4_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp

        isBilinear = False
        cquad_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp
        ctria3_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp
        ctria6_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp
        ctriar_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        if self.element_type == 74:
            msg = ctria3_msg
            nnodes = 3
            is_bilinear = False
        elif  self.element_type == 144:
            msg = cquad4_bilinear_msg
            nnodes = 4
            is_bilinear = True
        else:
            raise NotImplementedError(self.element_type)
        return msg, nnodes, is_bilinear


class RealPlateStrainArray(RealPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def isStress(self):
        return False
    def isStrain(self):
        return True

    def get_headers(self):
        if self.isFiberDistance():
            fd = 'fiber_distance'
        else:
            fd = 'fiber_curvature'

        if self.isVonMises():
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = [fd, 'exx', 'eyy', 'exy', 'angle', 'emax', 'emin', ovm]
        return headers

    def _get_msgs(self):
        if self.isVonMises():
            von_mises = 'VON MISES'
        else:
            von_mises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]

        #=============================

        isBilinear = False
        cquad4_msg = ['                         STRAINS   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp
        cquad8_msg = ['                         STRAINS   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        ## TODO: can cquad8s be bilinear???
        isBilinear = True
        cquadr_bilinear_msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
        cquad4_bilinear_msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp

        isBilinear = False
        cquad_msg = ['                         STRAINS   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp
        ctria3_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp
        ctria6_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp
        ctriar_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        is_bilinear = False
        if self.element_type == 74:
            msg = ctria3_msg
            nnodes = 3
        elif  self.element_type == 144:
            msg = cquad4_bilinear_msg
            nnodes = 4
            is_bilinear = True
        else:
            raise NotImplementedError(self.element_type)
        return msg, nnodes, is_bilinear

class RealPlateStress(StressObject):
    """
    ::

      ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)
        ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
            6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  4.227345E+03
                         1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  5.739119E+03


                           S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)          MAX
        ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR
            6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  4.227345E+03
                         1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  5.739119E+03
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
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.ovmShear = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                self.addNewNode = self.addNewNodeSort1
        else:
            raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2
            #self.addNewNode = self.addNewNodeSort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, fiberCurvature, oxx, oyy, txy, angle, '
                   'majorP, minorP, ovmShear\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType == 'CTRIA3':
                    (eType, eid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                                 f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'CEN/3': [f1, f2]}
                    self.oxx[eid] = {'CEN/3': [ox1, ox2]}
                    self.oyy[eid] = {'CEN/3': [oy1, oy2]}
                    self.txy[eid] = {'CEN/3': [txy1, txy2]}
                    self.angle[eid] = {'CEN/3': [angle1, angle2]}
                    self.majorP[eid] = {'CEN/3': [o11, o12]}
                    self.minorP[eid] = {'CEN/3': [o21, o22]}
                    self.ovmShear[eid] = {'CEN/3': [ovm1, ovm2]}
                elif eType == 'CQUAD4':
                    #assert len(line)==19,'len(line)=%s' %(len(line))
                    if len(line) == 19:  # Centroid - bilinear
                        (eType, eid, nid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                                          f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.oxx[eid] = {nid: [ox1, ox2]}
                        self.oyy[eid] = {nid: [oy1, oy2]}
                        self.txy[eid] = {nid: [txy1, txy2]}
                        self.angle[eid] = {nid: [angle1, angle2]}
                        self.majorP[eid] = {nid: [o11, o12]}
                        self.minorP[eid] = {nid: [o21, o22]}
                        self.ovmShear[eid] = {nid: [ovm1, ovm2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                                     f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        nid = 'CEN/4'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.oxx[eid] = {nid: [ox1, ox2]}
                        self.oyy[eid] = {nid: [oy1, oy2]}
                        self.txy[eid] = {nid: [txy1, txy2]}
                        self.angle[eid] = {nid: [angle1, angle2]}
                        self.majorP[eid] = {nid: [o11, o12]}
                        self.minorP[eid] = {nid: [o21, o22]}
                        self.ovmShear[eid] = {nid: [ovm1, ovm2]}
                    elif len(line) == 17:  # Bilinear
                        (nid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                              f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        self.fiberCurvature[eid][nid] = [f1, f2]
                        self.oxx[eid][nid] = [ox1, ox2]
                        self.oyy[eid][nid] = [oy1, oy2]
                        self.txy[eid][nid] = [txy1, txy2]
                        self.angle[eid][nid] = [angle1, angle2]
                        self.majorP[eid][nid] = [o11, o12]
                        self.minorP[eid][nid] = [o21, o22]
                        self.ovmShear[eid][nid] = [ovm1, ovm2]
                    else:
                        assert len(line) == 19, 'line=%s len(line)=%s' % (line, len(line))
                        raise NotImplementedError()
                else:
                    msg = 'eType=%r not supported...' % eType
                    raise NotImplementedError(msg)
            return

        dt = transient[1]
        if dt not in self.oxx:
            self.fiberCurvature[dt] = {}
            self.oxx[dt] = {}
            self.oyy[dt] = {}
            self.txy[dt] = {}
            self.angle[dt] = {}
            self.majorP[dt] = {}
            self.minorP[dt] = {}
            self.ovmShear[dt] = {}

        eType = data[0][0]
        for line in data:
            if eType == 'CTRIA3':
                (eType, eid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                             f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                self.eType[eid] = eType
                self.fiberCurvature[eid] = {'CEN/3': [f1, f2]}
                self.oxx[dt][eid] = {'CEN/3': [ox1, ox2]}
                self.oyy[dt][eid] = {'CEN/3': [oy1, oy2]}
                self.txy[dt][eid] = {'CEN/3': [txy1, txy2]}
                self.angle[dt][eid] = {'CEN/3': [angle1, angle2]}
                self.majorP[dt][eid] = {'CEN/3': [o11, o12]}
                self.minorP[dt][eid] = {'CEN/3': [o21, o22]}
                self.ovmShear[dt][eid] = {'CEN/3': [ovm1, ovm2]}
            elif eType == 'CQUAD4':
                if len(line) == 18:  # Centroid
                    (eType, eid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                                 f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                    nid = 'CEN/4'
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {nid: [f1, f2]}
                    self.oxx[dt][eid] = {nid: [ox1, ox2]}
                    self.oyy[dt][eid] = {nid: [oy1, oy2]}
                    self.txy[dt][eid] = {nid: [txy1, txy2]}
                    self.angle[dt][eid] = {nid: [angle1, angle2]}
                    self.majorP[dt][eid] = {nid: [o11, o12]}
                    self.minorP[dt][eid] = {nid: [o21, o22]}
                    self.ovmShear[dt][eid] = {nid: [ovm1, ovm2]}
                elif len(line) == 19:  # Centroid - bilinear
                    (eType, eid, nid, f1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                                      f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {nid: [f1, f2]}
                    self.oxx[dt][eid] = {nid: [ox1, ox2]}
                    self.oyy[dt][eid] = {nid: [oy1, oy2]}
                    self.txy[dt][eid] = {nid: [txy1, txy2]}
                    self.angle[dt][eid] = {nid: [angle1, angle2]}
                    self.majorP[dt][eid] = {nid: [o11, o12]}
                    self.minorP[dt][eid] = {nid: [o21, o22]}
                    self.ovmShear[dt][eid] = {nid: [ovm1, ovm2]}
                elif len(line) == 17:  # Bilinear node
                    (nid, f1, ox1, oy1, oxy1, angle1, o11, o21, ovm1,
                          f2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                    self.fiberCurvature[eid][nid] = [f1, f2]
                    self.oxx[dt][eid][nid] = [ox1, ox2]
                    self.oyy[dt][eid][nid] = [oy1, oy2]
                    self.txy[dt][eid][nid] = [oxy1, txy2]
                    self.angle[dt][eid][nid] = [angle1, angle2]
                    self.majorP[dt][eid][nid] = [o11, o12]
                    self.minorP[dt][eid][nid] = [o21, o22]
                    self.ovmShear[dt][eid][nid] = [ovm1, ovm2]
                else:
                    msg = 'line=%r not supported...len=%i' % (line, len(line))
                    raise NotImplementedError(msg)
            else:
                msg = 'eType=%r not supported...' % eType
                raise NotImplementedError(msg)

    def delete_transient(self, dt):
        #del self.fiberCurvature[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.txy[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.ovmShear[dt]

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
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.ovmShear[dt] = {}

    def add_new_eid(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #if eid in self.oxx:
            #return self.add(dt,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)

        #assert eid not in self.oxx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [fd]}
        assert isinstance(eid, int)
        self.oxx[eid] = {nodeID: [oxx]}
        self.oyy[eid] = {nodeID: [oyy]}
        self.txy[eid] = {nodeID: [txy]}
        self.angle[eid] = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.ovmShear[eid] = {nodeID: [ovm]}
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print msg
        if nodeID == 0:
            raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print('add_new_eid_sort1', msg)
        #print(self.oxx)
        #msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g major=%g vm=%g" % (dt, eid, nodeID, fd, oxx, majorP, ovm)
        #print msg
        #if eid in self.ovmShear[dt]:
        #    return self.add(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)

        if 0:
            if dt in self.oxx and eid in self.oxx[dt]:  # SOL200, erase the old result
                #nid = nodeID
                msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s major=%s vm=%s" %(dt,eid,nodeID,
                                                                             str(self.fiberCurvature[eid][nodeID]),
                                                                             str(self.oxx[dt][eid][nodeID]),
                                                                             str(self.majorP[dt][eid][nodeID]),
                                                                             str(self.ovmShear[dt][eid][nodeID]))
                #print("************", msg)
                self.delete_transient(dt)
                self.add_new_transient(dt)

        assert isinstance(eid, int)
        self.eType[eid] = eType
        if dt not in self.oxx:
            self.add_new_transient(dt)
        self.fiberCurvature[eid] = {nodeID: [fd]}
        self.oxx[dt][eid] = {nodeID: [oxx]}
        self.oyy[dt][eid] = {nodeID: [oyy]}
        self.txy[dt][eid] = {nodeID: [txy]}
        self.angle[dt][eid] = {nodeID: [angle]}
        self.majorP[dt][eid] = {nodeID: [majorP]}
        self.minorP[dt][eid] = {nodeID: [minorP]}
        self.ovmShear[dt][eid] = {nodeID: [ovm]}
        #print msg
        if nodeID == 0:
            msg = "dt=%s eid=%s nodeID=%s " % (dt, eid, nodeID) #fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(dt,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
            raise ValueError(msg)

    def add(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        assert isinstance(eid, int)
        self.fiberCurvature[eid][nodeID].append(fd)
        self.oxx[eid][nodeID].append(oxx)
        self.oyy[eid][nodeID].append(oyy)
        self.txy[eid][nodeID].append(txy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.ovmShear[eid][nodeID].append(ovm)
        if nodeID == 0:
            raise ValueError(msg)

    def add_sort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" % (dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print('add_sort1')
        #print msg
        #print self.oxx
        #print self.fiberCurvatrure
        assert eid is not None
        self.fiberCurvature[eid][nodeID].append(fd)
        self.oxx[dt][eid][nodeID].append(oxx)
        self.oyy[dt][eid][nodeID].append(oyy)
        self.txy[dt][eid][nodeID].append(txy)
        self.angle[dt][eid][nodeID].append(angle)
        self.majorP[dt][eid][nodeID].append(majorP)
        self.minorP[dt][eid][nodeID].append(minorP)
        self.ovmShear[dt][eid][nodeID].append(ovm)
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNode(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        assert eid is not None
        #assert nodeID not in self.oxx[eid]
        self.fiberCurvature[eid][nodeID] = [fd]
        self.oxx[eid][nodeID] = [oxx]
        self.oyy[eid][nodeID] = [oyy]
        self.txy[eid][nodeID] = [txy]
        self.angle[eid][nodeID] = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.ovmShear[eid][nodeID] = [ovm]
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNodeSort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #print self.oxx
        #print('addNewNodeSort1')
        assert eid is not None
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #print(msg)
        #assert nodeID not in self.oxx[dt][eid]
        self.fiberCurvature[eid][nodeID] = [fd]
        self.oxx[dt][eid][nodeID] = [oxx]
        self.oyy[dt][eid][nodeID] = [oyy]
        self.txy[dt][eid][nodeID] = [txy]
        self.angle[dt][eid][nodeID] = [angle]
        self.majorP[dt][eid][nodeID] = [majorP]
        self.minorP[dt][eid][nodeID] = [minorP]
        self.ovmShear[dt][eid][nodeID] = [ovm]
        if nodeID == 0:
            raise ValueError(msg)

    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']
        headers += ['oxx', 'oyy', 'txy', 'majorP', 'minorP']
        if self.isVonMises():
            headers.append('oVonMises')
        else:
            headers.append('maxShear')
        return headers

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f, is_mag_phase)

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]

        triMsg = None
        tri6Msg = None
        trirMsg = None
        quadMsg = None
        quad8Msg = None
        quadrMsg = None
        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[kkey].keys()
            isBilinear = True
            quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CQUAD8' in eTypes:
            quad8Msg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        if 'CQUADR' in eTypes:
            qkey = eTypes.index('CQUADR')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[kkey].keys()
            isBilinear = True
            quadrMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        for eType in typesOut:
            eids = orderedETypes[eType]
            if eids:
                eids.sort()
                #print "eType = ",eType
                #print "eids = ",eids
                #print "eType = ",eType
                msgPack = msgPacks[eType]

                msg += header + msgPack
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for eid in eids:
                            out = self._write_f06_quad4_bilinear(eid, 4, 'CEN/4')
                            msg.append(out)
                    else:
                        for eid in eids:
                            out = self._write_f06_tri3(eid)
                            msg.append(out)
                elif eType in ['CTRIA3']:
                    for eid in eids:
                        out = self._write_f06_tri3(eid)
                        msg.append(out)
                elif eType in ['CQUAD8']:
                    for eid in eids:
                        out = self._write_f06_quad4_bilinear(eid, 4, 'CEN/8')
                        msg.append(out)
                elif eType in ['CQUADR']:
                    for eid in eids:
                        out = self._write_f06_quad4_bilinear(eid, 4, 'CEN/4')
                        msg.append(out)
                elif eType in ['CTRIAR']:
                    for eid in eids:
                        out = self._write_f06_quad4_bilinear(eid, 3, 'CEN/3')
                        msg.append(out)
                elif eType in ['CTRIA6']:
                    for eid in eids:
                        out = self._write_f06_quad4_bilinear(eid, 3, 'CEN/6')
                        msg.append(out)
                else:
                    raise NotImplementedError('eType = %r' % eType)

                msg.append(pageStamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        return page_num - 1

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]

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
            #print(self.eType)
            kkey = self.eType.keys()[qkey]
            try:
                ekey = self.oxx[dt][kkey].keys()
            except KeyError:
                assert dt in self.oxx, 'dt=%r not in oxx' % dt
                assert kkey in self.oxx[dt], 'kkey=%r not in oxx[%r]' % (kkey, dt)
                raise
            isBilinear = True
            quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CQUAD8' in eTypes:
            quad8Msg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        if 'CQUADR' in eTypes:
            qkey = eTypes.index('CQUADR')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[dt][kkey].keys()
            isBilinear = True
            quadrMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        dts = self.oxx.keys()
        dts.sort()
        if isinstance(dts[0], int):
            dt_msg =  ' %s = %%-10i\n' % self.data_code['name']
        else:
            dt_msg =  ' %s = %%10.4E\n' % self.data_code['name']
        for eType in typesOut:
            #print "eType = ",eType
            eids = orderedETypes[eType]
            #print "eids = ",eids
            if eids:
                msgPack = msgPacks[eType]
                eids.sort()
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for dt in dts:
                            header[1] = dt_msg % dt
                            msg += header + msgPack
                            for eid in eids:
                                out = self._write_f06_quad4_bilinear_transient(dt, eid, 4, 'CEN/4')
                                msg.append(out)
                    else:
                        for dt in dts:
                            header[1] = dt_msg % dt
                            msg += header + msgPack
                            for eid in eids:
                                out = self._write_f06_tri3_transient(dt, eid)
                                msg.append(out)
                elif eType in ['CTRIA3']:
                    for dt in dts:
                        header[1] = dt_msg % dt
                        msg += header + msgPack
                        for eid in eids:
                            out = self._write_f06_tri3_transient(dt, eid)
                            msg.append(out)
                elif eType in ['CTRIAR']:
                    for dt in dts:
                        header[1] = dt_msg % dt
                        msg += header + msgPack
                        for eid in eids:
                            out = self._write_f06_quad4_bilinear_transient(dt, eid, 3, 'CEN/3')
                            msg.append(out)
                elif eType in ['CTRIA6']:
                    for dt in dts:
                        header[1] = dt_msg % dt
                        msg += header + msgPack
                        for eid in eids:
                            out = self._write_f06_quad4_bilinear_transient(dt, eid, 3, 'CEN/6')
                            msg.append(out)
                elif eType in ['CQUAD8']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                        msg += header + msgPack
                        for eid in eids:
                            out = self._write_f06_quad4_bilinear_transient(dt, eid, 5, 'CEN/8')
                            msg.append(out)
                elif eType in ['CQUADR']:
                    if isBilinear:
                        for dt in dts:
                            header[1] = dt_msg % dt
                            msg += header + msgPack
                            for eid in eids:
                                out = self._write_f06_quad4_bilinear_transient(dt, eid, 4, 'CEN/4')
                                msg.append(out)
                    else:
                        for dt in dts:
                            header[1] = dt_msg % dt
                            msg += header + msgPack
                            for eid in eids:
                                out = self._write_f06_tri3_transient(dt, eid)
                                msg.append(out)
                else:
                    raise NotImplementedError('eType = %r' % eType)  # CQUAD8, CTRIA6

                msg.append(pageStamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        return page_num - 1

    def _write_f06_quad4_bilinear(self, eid, n, cen):
        msg = ['']
        k = self.oxx[eid].keys()
        #cen = 'CEN/' + str(n)
        k.remove(cen)
        k.sort()
        for nid in [cen] + k:
            for iLayer in range(len(self.oxx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[eid][nid][iLayer]
                oyy = self.oyy[eid][nid][iLayer]
                txy = self.txy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                ovm = self.ovmShear[eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], is_all_zeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], is_all_zeros) = writeFloats8p4F([angle])

                if nid == cen and iLayer == 0:
                    msg.append('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % (eid, cen, fd, oxx, oyy, txy, angle, major, minor, ovm))
                elif iLayer == 0:
                    msg.append('   %8s %8i  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % ('', nid, fd, oxx, oyy, txy, angle, major, minor, ovm))
                elif iLayer == 1:
                    msg.append('   %8s %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n\n' % ('', '', fd, oxx, oyy, txy, angle, major, minor, ovm))
                else:
                    raise Exception('Invalid option for cquad4')
        return ''.join(msg)

    def _write_f06_quad4_bilinear_transient(self, dt, eid, n, cen):
        msg = ['']
        k = self.oxx[dt][eid].keys()
        #cen = 'CEN/' + str(n)
        k.remove(cen)
        k.sort()
        for nid in [cen] + k:
            for iLayer in range(len(self.oxx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[dt][eid][nid][iLayer]
                oyy = self.oyy[dt][eid][nid][iLayer]
                txy = self.txy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                ovm = self.ovmShear[dt][eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], is_all_zeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], is_all_zeros) = writeFloats8p4F([angle])

                if nid == cen and iLayer == 0:
                    msg.append('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % (eid, cen, fd, oxx, oyy, txy, angle, major, minor, ovm))
                elif iLayer == 0:
                    msg.append('   %8s %8i  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % ('', nid, fd, oxx, oyy, txy, angle, major, minor, ovm))
                elif iLayer == 1:
                    msg.append('   %8s %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n\n' % ('', '', fd, oxx, oyy, txy, angle, major, minor, ovm))
                else:
                    raise RuntimeError('Invalid option for cquad4')
        return ''.join(msg)

    def _write_f06_tri3(self, eid):
        msg = ['']
        oxxNodes = self.oxx[eid].keys()
        for nid in sorted(oxxNodes):
            for iLayer in range(len(self.oxx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[eid][nid][iLayer]
                oyy = self.oyy[eid][nid][iLayer]
                txy = self.txy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                ovm = self.ovmShear[eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], is_all_zeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], is_all_zeros) = writeFloats8p4F([angle])

                if iLayer == 0:
                    msg.append('0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % (eid, fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip()))
                else:
                    msg.append('   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % ('', fd, oxx, oyy, txy, angle, major, minor, ovm.rstrip()))
        return ''.join(msg)

    def _write_f06_tri3_transient(self, dt, eid):
        msg = ['']
        oxxNodes = self.oxx[dt][eid].keys()
        for nid in sorted(oxxNodes):
            for iLayer in range(len(self.oxx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                oxx = self.oxx[dt][eid][nid][iLayer]
                oyy = self.oyy[dt][eid][nid][iLayer]
                txy = self.txy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                ovm = self.ovmShear[dt][eid][nid][iLayer]
                ([fd, oxx, oyy, txy, major, minor, ovm], is_all_zeros) = writeFloats13E([fd, oxx, oyy, txy, major, minor, ovm])
                ([angle], is_all_zeros) = writeFloats8p4F([angle])

                if iLayer == 0:
                    msg.append('0  %6i   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %s\n' % (eid, fd, oxx, oyy, txy, angle, major, minor, ovm))
                else:
                    msg.append('   %6s   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %s\n' % ('', fd, oxx, oyy, txy, angle, major, minor, ovm))
        return ''.join(msg)


class RealPlateStrain(StrainObject):
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

        #print self.data_code
        self.fiberCurvature = {}
        self.exx = {}
        self.eyy = {}
        self.exy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.evmShear = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                self.addNewNode = self.addNewNodeSort1
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        else:
            raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

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
        msg.append('  eType, fiberCurvature, exx, eyy, exy, angle, '
                   'majorP, minorP, evmShear\n')
        return msg

    def add_f06_data(self, data, transient, data_code=None):
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType == 'CTRIA3':
                    (eType, eid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                                 f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'CEN/3': [f1, f2]}
                    self.exx[eid] = {'CEN/3': [ex1, ex2]}
                    self.eyy[eid] = {'CEN/3': [ey1, ey2]}
                    self.exy[eid] = {'CEN/3': [exy1, exy2]}
                    self.angle[eid] = {'CEN/3': [angle1, angle2]}
                    self.majorP[eid] = {'CEN/3': [e11, e12]}
                    self.minorP[eid] = {'CEN/3': [e21, e22]}
                    self.evmShear[eid] = {'CEN/3': [evm1, evm2]}
                elif eType == 'CQUAD4':
                    #assert len(line)==19,'len(line)=%s' %(len(line))
                    #print line
                    if len(line) == 19:  # Centroid - bilinear
                        (
                            eType, eid, nid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                                             f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.exx[eid] = {nid: [ex1, ex2]}
                        self.eyy[eid] = {nid: [ey1, ey2]}
                        self.exy[eid] = {nid: [exy1, exy2]}
                        self.angle[eid] = {nid: [angle1, angle2]}
                        self.majorP[eid] = {nid: [e11, e12]}
                        self.minorP[eid] = {nid: [e21, e22]}
                        self.evmShear[eid] = {nid: [evm1, evm2]}
                    elif len(line) == 18:  # Centroid
                        (eType, eid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                                     f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        nid = 'CEN/4'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid: [f1, f2]}
                        self.exx[eid] = {nid: [ex1, ex2]}
                        self.eyy[eid] = {nid: [ey1, ey2]}
                        self.exy[eid] = {nid: [exy1, exy2]}
                        self.angle[eid] = {nid: [angle1, angle2]}
                        self.majorP[eid] = {nid: [e11, e12]}
                        self.minorP[eid] = {nid: [e21, e22]}
                        self.evmShear[eid] = {nid: [evm1, evm2]}
                    elif len(line) == 17:  # Bilinear node
                        (nid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                              f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        self.fiberCurvature[eid][nid] = [f1, f2]
                        self.exx[eid][nid] = [ex1, ex2]
                        self.eyy[eid][nid] = [ey1, ey2]
                        self.exy[eid][nid] = [exy1, exy2]
                        self.angle[eid][nid] = [angle1, angle2]
                        self.majorP[eid][nid] = [e11, e12]
                        self.minorP[eid][nid] = [e21, e22]
                        self.evmShear[eid][nid] = [evm1, evm2]
                    else:
                        #assert len(line) == 19, 'len(line)=%s' % len(line)
                        msg = 'line=%r not supported...len=%i' % (line, len(line))
                        raise NotImplementedError(msg)
                else:
                    msg = 'line=%s not supported...' % line
                    raise NotImplementedError(msg)
            return
        eType = data[0][0]
        assert 'name' in self.data_code, self.data_code

        dt = transient[1]
        if dt not in self.exx:
            self.exx[dt] = {}
            self.eyy[dt] = {}
            self.exy[dt] = {}
            self.angle[dt] = {}
            self.majorP[dt] = {}
            self.minorP[dt] = {}
            self.evmShear[dt] = {}

        for line in data:
            eType = data[0][0]
            if eType == 'CTRIA3':
                (eType, eid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                             f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                self.eType[eid] = eType
                self.fiberCurvature[eid] = {'CEN/3': [f1, f2]}
                self.exx[dt][eid] = {'CEN/3': [ex1, ex2]}
                self.eyy[dt][eid] = {'CEN/3': [ey1, ey2]}
                self.exy[dt][eid] = {'CEN/3': [exy1, exy2]}
                self.angle[dt][eid] = {'CEN/3': [angle1, angle2]}
                self.majorP[dt][eid] = {'CEN/3': [e11, e12]}
                self.minorP[dt][eid] = {'CEN/3': [e21, e22]}
                self.evmShear[dt][eid] = {'CEN/3': [evm1, evm2]}
            elif eType == 'CQUAD4':
                if len(line) == 19:  # Centroid - bilinear
                    (
                        eType, eid, nid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                                         f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {nid: [f1, f2]}
                    self.exx[dt][eid] = {nid: [ex1, ex2]}
                    self.eyy[dt][eid] = {nid: [ey1, ey2]}
                    self.exy[dt][eid] = {nid: [exy1, exy2]}
                    self.angle[dt][eid] = {nid: [angle1, angle2]}
                    self.majorP[dt][eid] = {nid: [e11, e12]}
                    self.minorP[dt][eid] = {nid: [e21, e22]}
                    self.evmShear[dt][eid] = {nid: [evm1, evm2]}
                elif len(line) == 18:  # Centroid
                    (
                        eType, eid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                        f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {nid: [f1, f2]}
                    self.exx[dt][eid] = {nid: [ex1, ex2]}
                    self.eyy[dt][eid] = {nid: [ey1, ey2]}
                    self.exy[dt][eid] = {nid: [exy1, exy2]}
                    self.angle[dt][eid] = {nid: [angle1, angle2]}
                    self.majorP[dt][eid] = {nid: [e11, e12]}
                    self.minorP[dt][eid] = {nid: [e21, e22]}
                    self.evmShear[dt][eid] = {nid: [evm1, evm2]}
                elif len(line) == 17:  # Bilinear node
                    (nid, f1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                          f2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    self.fiberCurvature[eid][nid] = [f1, f2]
                    self.exx[dt][eid][nid] = [ex1, ex2]
                    self.eyy[dt][eid][nid] = [ey1, ey2]
                    self.exy[dt][eid][nid] = [exy1, exy2]
                    self.angle[dt][eid][nid] = [angle1, angle2]
                    self.majorP[dt][eid][nid] = [e11, e12]
                    self.minorP[dt][eid][nid] = [e21, e22]
                    self.evmShear[dt][eid][nid] = [evm1, evm2]
                else:
                    msg = 'line=%r not supported...len=%i' % (line, len(line))
                    raise NotImplementedError(msg)
            else:
                msg = 'eType=%r is not supported...' % eType
                raise NotImplementedError(msg)

    def delete_transient(self, dt):
        #del self.fiberCurvature[dt]
        del self.exx[dt]
        del self.eyy[dt]
        del self.exy[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.evmShear[dt]

    def get_transients(self):
        k = self.exx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #self.fiberCurvature[dt] = {}
        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.exy[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.evmShear[dt] = {}

    def add_new_eid(self, eType, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)

        #if nodeID != 'C':  # centroid
            #assert 0 < nodeID < 1000000000, 'nodeID=%s %s' % (nodeID, msg)

        #if eid in self.exx:
            #return self.add(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #assert eid not in self.exx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[eid] = {nodeID: [exx]}
        self.eyy[eid] = {nodeID: [eyy]}
        self.exy[eid] = {nodeID: [exy]}
        self.angle[eid] = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.evmShear[eid] = {nodeID: [evm]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add_new_eid_sort1(self, eType, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #print msg

        #if nodeID != 'C':  # centroid
            #assert 0 < nodeID < 1000000000, 'nodeID=%s %s' % (nodeID, msg)

        if dt not in self.exx:
            self.add_new_transient(dt)
        #if eid in self.evmShear[dt]:  # SOL200, erase the old result
            #nid = nodeID
            #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s major=%s vm=%s" %(dt,eid,nodeID,str(self.fiberCurvature[eid][nid]),str(self.oxx[dt][eid][nid]),str(self.majorP[dt][eid][nid]),str(self.ovmShear[dt][eid][nid]))
            #self.delete_transient(dt)
            #self.add_new_transient()

        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[dt][eid] = {nodeID: [exx]}
        self.eyy[dt][eid] = {nodeID: [eyy]}
        self.exy[dt][eid] = {nodeID: [exy]}
        self.angle[dt][eid] = {nodeID: [angle]}
        self.majorP[dt][eid] = {nodeID: [majorP]}
        self.minorP[dt][eid] = {nodeID: [minorP]}
        self.evmShear[dt][eid] = {nodeID: [evm]}
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add(self, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        #if nodeID != 'C':  # centroid
            #assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)
        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[eid][nodeID].append(exx)
        self.eyy[eid][nodeID].append(eyy)
        self.exy[eid][nodeID].append(exy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.evmShear[eid][nodeID].append(evm)
        if nodeID == 0:
            raise ValueError(msg)

    def add_sort1(self, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        #if nodeID != 'C':  # centroid
            #assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)

        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[dt][eid][nodeID].append(exx)
        self.eyy[dt][eid][nodeID].append(eyy)
        self.exy[dt][eid][nodeID].append(exy)
        self.angle[dt][eid][nodeID].append(angle)
        self.majorP[dt][eid][nodeID].append(majorP)
        self.minorP[dt][eid][nodeID].append(minorP)
        self.evmShear[dt][eid][nodeID].append(evm)
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNode(self, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print self.oxx
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        assert nodeID not in self.exx[eid], msg
        self.fiberCurvature[eid][nodeID] = [curvature]
        self.exx[eid][nodeID] = [exx]
        self.eyy[eid][nodeID] = [eyy]
        self.exy[eid][nodeID] = [exy]
        self.angle[eid][nodeID] = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.evmShear[eid][nodeID] = [evm]
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def addNewNodeSort1(self, dt, eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #print self.oxx
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, nodeID, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #assert nodeID not in self.exx[eid], msg
        self.fiberCurvature[eid][nodeID] = [curvature]
        self.exx[dt][eid][nodeID] = [exx]
        self.eyy[dt][eid][nodeID] = [eyy]
        self.exy[dt][eid][nodeID] = [exy]
        self.angle[dt][eid][nodeID] = [angle]
        self.majorP[dt][eid][nodeID] = [majorP]
        self.minorP[dt][eid][nodeID] = [minorP]
        self.evmShear[dt][eid][nodeID] = [evm]
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']

        headers += ['exx', 'eyy', 'exy', 'eMajor', 'eMinor']
        if self.isVonMises():
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f)

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER                STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % (vonMises)]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % (vonMises)]

        quadMsg = None
        quad8Msg = None
        quadrMsg = None
        triMsg = None
        tri6Msg = None
        trirMsg = None

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            ekey = self.exx[kkey].keys()
            isBilinear = True
            quadMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CQUAD8' in eTypes:
            qkey = eTypes.index('CQUAD8')
            kkey = self.eType.keys()[qkey]
            ekey = self.exx[kkey].keys()
            isBilinear = True
            quad8Msg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quad8Msg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + triMsgTemp

        if 'CQUADR' in eTypes:
            qkey = eTypes.index('CQUADR')
            kkey = self.eType.keys()[qkey]
            ekey = self.exx[kkey].keys()
            isBilinear = True
            quadrMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadrMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        msg = []
        for eType in typesOut:
            eids = orderedETypes[eType]
            if eids:
                #print "eids = ",eids
                #print "eType = ",eType
                msgPack = msgPacks[eType]
                eids.sort()
                msg += header + msgPack
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for eid in eids:
                            out = self.writeF06_Quad4_Bilinear(eid, 4, 'CEN/4')
                            msg.append(out)
                    else:
                        for eid in eids:
                            out = self.writeF06_Tri3(eid)
                            msg.append(out)
                elif eType in ['CTRIA3']:
                    for eid in eids:
                        out = self.writeF06_Tri3(eid)
                        msg.append(out)
                elif eType in ['CQUAD8']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 4, 'CEN/8')
                        msg.append(out)
                elif eType in ['CQUADR']:
                    if isBilinear:
                        for eid in eids:
                            out = self.writeF06_Quad4_Bilinear(eid, 4, 'CEN/4')
                            msg.append(out)
                    else:
                        for eid in eids:
                            out = self.writeF06_Tri3(eid)
                            msg.append(out)
                elif eType in ['CTRIAR']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 3, 'CEN/3')
                        msg.append(out)
                elif eType in ['CTRIA6']:
                    for eid in eids:
                        out = self.writeF06_Quad4_Bilinear(eid, 3, 'CEN/6')
                        msg.append(out)
                else:
                    raise NotImplementedError('eType = %r' % eType)
                msg.append(pageStamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
        return page_num - 1

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER                STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % vonMises]
            triMsgTemp = ['  ELEMENT      FIBER               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % vonMises]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % vonMises]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % vonMises]

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
            isBilinear = True
            quadMsg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if nLayers == 1:
                isBilinear = False
                quadMsg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + triMsgTemp

        if 'CQUAD8' in eTypes:
            quad8Msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + quadMsgTemp

        if 'CQUADR' in eTypes:
            qkey = eTypes.index('CQUADR')
            kkey = self.eType.keys()[qkey]
            dt = self.exx.keys()[0]
            ekey = self.exx[dt][kkey].keys()
            isBilinear = True
            quadrMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quadMsgTemp
            if len(ekey) == 1:
                isBilinear = False
                quadrMsg = header + ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + triMsgTemp

        if 'CTRIA3' in eTypes:
            triMsg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + triMsgTemp

        if 'CTRIA6' in eTypes:
            tri6Msg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + triMsgTemp

        if 'CTRIAR' in eTypes:
            trirMsg = header + ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + triMsgTemp

        msg = []
        msgPacks = {'CTRIA3': triMsg,
                    'CTRIA6': tri6Msg,
                    'CTRIAR': trirMsg,
                    'CQUAD4': quadMsg,
                    'CQUAD8': quad8Msg,
                    'CQUADR': quadrMsg, }

        validTypes = ['CTRIA3', 'CTRIA6', 'CTRIAR', 'CQUAD4',
                      'CQUAD8', 'CQUADR']
        (typesOut, orderedETypes) = self.getOrderedETypes(validTypes)

        dts = self.exx.keys()
        dts.sort()
        if isinstance(dts[0], int):
            dt_msg =  ' %s = %%-10i\n' % self.data_code['name']
        else:
            dt_msg =  ' %s = %%10.4E\n' % self.data_code['name']

        for eType in typesOut:
            eids = orderedETypes[eType]
            if eids:
                msg = []
                msg_pack = msgPacks[eType]
                eids.sort()
                if eType in ['CQUAD4']:
                    if isBilinear:
                        for dt in dts:
                            header[1] = dt_msg % dt
                            msg.append('\n'.join(header + msg_pack))
                            for eid in eids:
                                out = self._write_f06_quad4_bilinear_transient(dt, eid, 4, 'CEN/4')
                                msg.append(out)
                            msg.append(pageStamp % page_num)
                            page_num += 1
                    else:
                        for dt in dts:
                            header[1] = dt_msg % dt
                            msg.append('\n'.join(header + msg_pack))
                            for eid in eids:
                                out = self._write_f06_tri3_transient(dt, eid)
                                msg.append(out)
                            msg.append(pageStamp % page_num)
                            page_num += 1
                elif eType in ['CTRIA3']:
                    for dt in dts:
                        header[1] = dt_msg % dt
                        msg.append('\n'.join(header + msg_pack))
                        for eid in eids:
                            out = self._write_f06_tri3_transient(dt, eid)
                            msg.append(out)
                        msg.append(pageStamp % page_num)
                        page_num += 1
                elif eType in ['CQUAD8']:
                    for dt in dts:
                        header[1] = dt_msg % dt
                        msg.append('\n'.join(header + msg_pack))
                        for eid in eids:
                            out = self._write_f06_quad4_bilinear_transient(dt, eid, 5, 'CEN/8')
                            msg.append(out)
                        msg.append(pageStamp % page_num)
                        page_num += 1
                elif eType in ['CQUADR']:
                    if isBilinear:
                        for dt in dts:
                            header[1] = dt_msg % dt
                            msg += header + msg_pack
                            for eid in eids:
                                out = self._write_f06_quad4_bilinear_transient(dt, eid, 4, 'CEN/4')
                                msg.append(out)
                    else:
                        for dt in dts:
                            header[1] = dt_msg % dt
                            msg += header + msg_pack
                            for eid in eids:
                                out = self._write_f06_tri3_transient(dt, eid)
                                msg.append(out)
                elif eType in ['CTRIAR']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                        msg.append('\n'.join(header + msg_pack))
                        for eid in eids:
                            out = self._write_f06_quad4_bilinear_transient(dt, eid, 3, 'CEN/3')
                            msg.append(out)
                        msg.append(pageStamp % page_num)
                        page_num += 1
                elif eType in ['CTRIA6']:
                    for dt in dts:
                        header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                        msg.append('\n'.join(header + msg_pack))
                        for eid in eids:
                            out = self._write_f06_quad4_bilinear_transient(dt, eid, 3, 'CEN/6')
                            msg.append(out)
                        msg.append(pageStamp % page_num)
                        page_num += 1
                else:
                    raise NotImplementedError('eType = %r' % eType)  # CQUAD8, CTRIA6
                f.write(''.join(msg))
        return page_num - 1

    def writeF06_Quad4_Bilinear(self, eid, n, cen):
        msg = ['']
        k = self.exx[eid].keys()
        #cen = 'CEN/' + str(n)
        k.remove(cen)
        k.sort()
        for nid in [cen] + k:
            for iLayer in range(len(self.exx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[eid][nid][iLayer]
                eyy = self.eyy[eid][nid][iLayer]
                exy = self.exy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                evm = self.evmShear[eid][nid][iLayer]
                ([fd, exx, eyy, exy, major, minor, evm], is_all_zeros) = writeFloats13E([fd, exx, eyy, exy, major, minor, evm])
                ([angle], is_all_zeros) = writeFloats8p4F([angle])

                if nid == cen and iLayer == 0:
                    msg.append('0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %s\n' % (eid, cen, fd, exx, eyy, exy, angle, major, minor, evm))
                elif iLayer == 0:
                    msg.append('   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %s\n' % ('', nid, fd, exx, eyy, exy, angle, major, minor, evm))
                elif iLayer == 1:
                    msg.append('   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %s\n\n' % ('', '', fd, exx, eyy, exy, angle, major, minor, evm))
                else:
                    raise RuntimeError('Invalid option for cquad4')
        return ''.join(msg)

    def _write_f06_quad4_bilinear_transient(self, dt, eid, n, cen):
        msg = ['']
        k = self.exx[dt][eid].keys()
        #cen = 'CEN/' + str(n)
        k.remove(cen)
        k.sort()
        for nid in [cen] + k:
            for iLayer in range(len(self.exx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[dt][eid][nid][iLayer]
                eyy = self.eyy[dt][eid][nid][iLayer]
                exy = self.exy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                evm = self.evmShear[dt][eid][nid][iLayer]

                ([fd, exx, eyy, exy, major, minor, evm], is_all_zeros) = writeFloats13E([fd, exx, eyy, exy, major, minor, evm])
                ([angle], is_all_zeros) = writeFloats8p4F([angle])

                if nid == cen and iLayer == 0:
                    msg.append('0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % (eid, cen, fd, exx, eyy, exy, angle, major, minor, evm.rstrip()))
                elif iLayer == 0:
                    msg.append('   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' % ('', nid, fd, exx, eyy, exy, angle, major, minor, evm.rstrip()))
                elif iLayer == 1:
                    msg.append('   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' % ('', '', fd, exx, eyy, exy, angle, major, minor, evm.rstrip()))
                else:
                    raise RuntimeError('Invalid option for cquad4')
        return ''.join(msg)

    def writeF06_Tri3(self, eid):
        msg = ['']
        k = self.exx[eid].keys()
        for nid in sorted(k):
            for iLayer in range(len(self.exx[eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[eid][nid][iLayer]
                eyy = self.eyy[eid][nid][iLayer]
                exy = self.exy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                evm = self.evmShear[eid][nid][iLayer]

                ([fd, exx, eyy, exy, major, minor, evm], is_all_zeros) = writeFloats13E([fd, exx, eyy, exy, major, minor, evm])
                ([angle], is_all_zeros) = writeFloats8p4F([angle])
                if iLayer == 0:
                    msg.append('0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % (eid, fd, exx, eyy, exy, angle, major, minor, evm.rstrip()))
                else:
                    msg.append('   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' % ('', fd, exx, eyy, exy, angle, major, minor, evm.rstrip()))
        return ''.join(msg)

    def _write_f06_tri3_transient(self, dt, eid):
        msg = ['']
        exxNodes = self.exx[dt][eid]
        #k = exxNodes.keys()
        for nid in sorted(exxNodes):
            for iLayer in range(len(self.exx[dt][eid][nid])):
                fd = self.fiberCurvature[eid][nid][iLayer]
                exx = self.exx[dt][eid][nid][iLayer]
                eyy = self.eyy[dt][eid][nid][iLayer]
                exy = self.exy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                evm = self.evmShear[dt][eid][nid][iLayer]

                ([fd, exx, eyy, exy, major, minor, evm], is_all_zeros) = writeFloats13E([fd, exx, eyy, exy, major, minor, evm])
                ([angle], is_all_zeros) = writeFloats8p4F([angle])
                if iLayer == 0:
                    msg.append('0  %6i   %13s     %13s  %13s  %13s   %8s   '
                               '%13s   %13s  %-s\n' % (eid, fd, exx, eyy, exy,
                                                    angle, major, minor, evm))
                else:
                    msg.append('   %6s   %13s     %13s  %13s  %13s   %8s   '
                               '%13s   %13s  %-s\n' % ('', fd, exx, eyy, exy,
                                                    angle, major, minor, evm))
        return ''.join(msg)
