#pylint disable=C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip, range
from itertools import count
from numpy import zeros, searchsorted, ravel, array_equal

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, writeFloats8p4F, _eigenvalue_header, get_key0

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

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def is_bilinear(self):
        if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
            return False
        elif self.element_type in [144, 64, 82, 70, 75]:  # CQUAD4
            return True
        else:
            raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))

    def build(self):
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type in [33, 74]:
            nnodes_per_element = 1
        elif self.element_type == 144:
            nnodes_per_element = 5
        elif self.element_type == 64:  # CQUAD8
            nnodes_per_element = 5
        elif self.element_type == 82:  # CQUADR
            nnodes_per_element = 5
        elif self.element_type == 70:  # CTRIAR
            nnodes_per_element = 4
        elif self.element_type == 75:  # CTRIA6
            nnodes_per_element = 4
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

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

        #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'element_node shape=%s table.shape=%s' % (self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += '(Eid, Nid)\n'
            for (eid, nid), (eid2, nid2) in zip(self.element_node, table.element_node):
                msg += '(%s, %s)    (%s, %s)\n' % (eid, nid, eid2, nid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm) = t1
                    (fiber_dist2, oxx2, oyy2, txy2, angle2, majorP2, minorP2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm,
                            fiber_dist2, oxx2, oyy2, txy2, angle2, majorP2, minorP2, ovm2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_new_eid(self, etype, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        self.add_new_eid_sort1(etype, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    def add_new_eid_sort1(self, etype, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        assert isinstance(eid, int), eid
        assert isinstance(node_id, int), node_id
        self._times[self.itime] = dt
        #assert self.itotal == 0, oxx
        self.element_node[self.itotal, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        self.itotal += 1
        self.ielement += 1

    def addNewNode(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        assert isinstance(node_id, int), node_id
        self.add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    def add(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        assert isinstance(node_id, int), node_id
        self.add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    def addNewNodeSort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        self.add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    def add_sort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        assert eid is not None, eid
        assert isinstance(node_id, int), node_id
        self.element_node[self.itotal, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
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
        nnodes = self.nnodes
        ntotal = self.ntotal
        nlayers = 2
        nelements = self.ntotal // self.nnodes // 2

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msgi = '  type=%s ntimes=%i nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, nlayers, ntotal)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, nelements, nnodes, nlayers, ntotal)
            ntimes_word = 1
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 str(', '.join(headers))))
        msg.append('  data.shape=%s\n' % str(self.data.shape))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg, nnodes, cen = _get_plate_msg(self)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        #cen_word = 'CEN/%i' % nnodes
        cen_word = cen
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
            fiber_dist = self.data[itime, :, 0]
            oxx = self.data[itime, :, 1]
            oyy = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            angle = self.data[itime, :, 4]
            majorP = self.data[itime, :, 5]
            minorP = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]

            # loop over all the elements
            for (i, eid, nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi) in zip(
                 count(), eids, nids, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
                ([fdi, oxxi, oyyi, txyi, major, minor, ovmi],
                 is_all_zeros) = writeFloats13E([fdi, oxxi, oyyi, txyi, major, minor, ovmi])
                ilayer = i % 2
                # tria3
                if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
                    if ilayer == 0:
                        f.write('0  %6i   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (eid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    else:
                        f.write('   %6s   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % ('', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))

                elif self.element_type in [64, 70, 75, 82, 144]:  # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4
                    # bilinear
                    if nid == 0 and ilayer == 0:  # CEN
                        f.write('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (eid, cen_word, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 0:
                        f.write('   %8s %8i  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % ('', nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 1:
                        f.write('   %8s %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n\n' % ('', '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                else:
                    raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealPlateStressArray(RealPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        if self.is_fiber_distance():
            fiber_dist = 'fiber_distance'
        else:
            fiber_dist = 'fiber_curvature'

        if self.is_von_mises():
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = [fiber_dist, 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', ovm]
        return headers


class RealPlateStrainArray(RealPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        if self.is_fiber_distance():
            fiber_dist = 'fiber_distance'
        else:
            fiber_dist = 'fiber_curvature'

        if self.is_von_mises():
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = [fiber_dist, 'exx', 'eyy', 'exy', 'angle', 'emax', 'emin', ovm]
        return headers


def _get_plate_msg(self):
    if self.is_von_mises():
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    if self.is_stress():
        if self.is_fiber_distance():
            quad_msg_temp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                            '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        else:
            quad_msg_temp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                             '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                            '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]

        #is_bilinear = False
        cquad4_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + tri_msg_temp
        cquad8_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + tri_msg_temp
        cquadr_msg = ['                        S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp

        ## TODO: can cquad8s be bilinear???
        #is_bilinear = True
        #cquadr_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quad_msg_temp
        cquad4_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quad_msg_temp

        #is_bilinear = False
        ctria3_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + tri_msg_temp
        ctria6_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + tri_msg_temp
        ctriar_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + tri_msg_temp
    else:
        if self.is_fiber_distance():
            quad_msg_temp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                            '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        else:
            quad_msg_temp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                            '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]

        #is_bilinear = False
        cquad4_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + tri_msg_temp
        cquad8_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + tri_msg_temp
        cquadr_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp

        ## TODO: can cquad8s be bilinear???
        #cquadr_bilinear_msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quad_msg_temp
        cquad4_bilinear_msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quad_msg_temp

        cquadr_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp
        ctria3_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + tri_msg_temp
        ctria6_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + tri_msg_temp
        ctriar_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + tri_msg_temp


    #is_bilinear = False
    if self.element_type == 74:
        msg = ctria3_msg
        nnodes = 3
        cen = 'CEN/3'
    elif self.element_type == 33:
        msg = cquad4_msg
        nnodes = 4
        cen = 'CEN/4'
    elif self.element_type == 144:
        msg = cquad4_bilinear_msg
        nnodes = 4
        #is_bilinear = True
        cen = 'CEN/4'
    elif self.element_type == 82:  # CQUADR
        msg = cquadr_msg
        nnodes = 4
        #is_bilinear = True
        cen = 'CEN/4'
    elif self.element_type == 64:  # CQUAD8
        msg = cquad8_msg
        nnodes = 4
        #is_bilinear = True
        cen = 'CEN/8'
    elif self.element_type == 75:  # CTRIA6
        msg = ctria6_msg
        nnodes = 3
        #is_bilinear = True
        cen = 'CEN/6'
    elif self.element_type == 70:  # CTRIAR
        msg = ctriar_msg
        nnodes = 3
        #is_bilinear = True
        cen = 'CEN/3'
    else:
        raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))
    return msg, nnodes, cen

#class RealPlateStress(StressObject):
    #"""
    #::

      #ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)
        #ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
            #6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  4.227345E+03
                         #1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  5.739119E+03


                           #S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      #ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)          MAX
        #ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR
            #6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  4.227345E+03
                         #1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  5.739119E+03
    #"""
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StressObject.__init__(self, data_code, isubcase)
        #self.eType = {}

        ##self.append_data_member('sCodes','s_code')
        ##print "self.s_code = ",self.s_code
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.fiberCurvature = {}
        #self.oxx = {}
        #self.oyy = {}
        #self.txy = {}
        #self.angle = {}
        #self.majorP = {}
        #self.minorP = {}
        #self.ovmShear = {}

        #self.dt = dt
        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
                #self.addNewNode = self.addNewNodeSort1
        #else:
            #raise NotImplementedError('SORT2')

    #def get_stats(self):
        #nelements = len(self.eType)

        #msg = self.get_data_code()
        #if self.nonlinear_factor is not None:  # transient
            #ntimes = len(self.oxx)
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
        #msg.append('  eType, fiberCurvature, oxx, oyy, txy, angle, '
                   #'majorP, minorP, ovmShear\n')
        #msg.append('  scode=%s stress_bits=%s is_stress=%s dist=%s vm=%s\n' % (
            #self.s_code, self.stress_bits, self.is_stress(), self.is_fiber_distance(), self.is_von_mises()))
        #return msg

    #def add_f06_data(self, data, transient):
        #if transient is None:
            #etype = data[0][0]
            #for line in data:
                #if etype == 'CTRIA3':
                    #cen = 0 # 'CEN/3'
                    #(etype, eid,
                     #fiber_dist1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                     #fiber_dist2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                    #self.eType[eid] = etype
                    #self.fiberCurvature[eid] = {cen : [fiber_dist1, fiber_dist2]}
                    #self.oxx[eid] = {cen : [ox1, ox2]}
                    #self.oyy[eid] = {cen : [oy1, oy2]}
                    #self.txy[eid] = {cen : [txy1, txy2]}
                    #self.angle[eid] = {cen : [angle1, angle2]}
                    #self.majorP[eid] = {cen : [o11, o12]}
                    #self.minorP[eid] = {cen : [o21, o22]}
                    #self.ovmShear[eid] = {cen : [ovm1, ovm2]}
                #elif etype == 'CQUAD4':
                    #if len(line) == 19:  # Centroid - bilinear
                        #(etype, eid, nid,
                         #fiber_dist1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                         #fiber_dist2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        #if nid == 'CEN/4':
                            #nid = 0
                        #assert isinstance(nid, int), nid
                        #self.eType[eid] = etype
                        #self.fiberCurvature[eid] = {nid : [fiber_dist1, fiber_dist2]}
                        #self.oxx[eid] = {nid : [ox1, ox2]}
                        #self.oyy[eid] = {nid : [oy1, oy2]}
                        #self.txy[eid] = {nid : [txy1, txy2]}
                        #self.angle[eid] = {nid : [angle1, angle2]}
                        #self.majorP[eid] = {nid : [o11, o12]}
                        #self.minorP[eid] = {nid : [o21, o22]}
                        #self.ovmShear[eid] = {nid : [ovm1, ovm2]}
                    #elif len(line) == 18:  # Centroid
                        #(etype, eid,
                         #fiber_dist1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                         #fiber_dist2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        #nid = 0 # CEN/4
                        #self.eType[eid] = etype
                        #self.fiberCurvature[eid] = {nid : [fiber_dist1, fiber_dist2]}
                        #self.oxx[eid] = {nid : [ox1, ox2]}
                        #self.oyy[eid] = {nid : [oy1, oy2]}
                        #self.txy[eid] = {nid : [txy1, txy2]}
                        #self.angle[eid] = {nid : [angle1, angle2]}
                        #self.majorP[eid] = {nid : [o11, o12]}
                        #self.minorP[eid] = {nid : [o21, o22]}
                        #self.ovmShear[eid] = {nid : [ovm1, ovm2]}
                    #elif len(line) == 17:  # Bilinear
                        #(nid,
                         #fiber_dist1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                         #fiber_dist2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                        #assert isinstance(nid, int), nid
                        #self.fiberCurvature[eid][nid] = [fiber_dist1, fiber_dist2]
                        #self.oxx[eid][nid] = [ox1, ox2]
                        #self.oyy[eid][nid] = [oy1, oy2]
                        #self.txy[eid][nid] = [txy1, txy2]
                        #self.angle[eid][nid] = [angle1, angle2]
                        #self.majorP[eid][nid] = [o11, o12]
                        #self.minorP[eid][nid] = [o21, o22]
                        #self.ovmShear[eid][nid] = [ovm1, ovm2]
                    #else:
                        #assert len(line) == 19, 'line=%s len(line)=%s' % (line, len(line))
                        #raise NotImplementedError()
                #else:
                    #msg = 'eType=%r not supported...' % etype
                    #raise NotImplementedError(msg)
            #return

        #dt = transient[1]
        #if dt not in self.oxx:
            #self.fiberCurvature[dt] = {}
            #self.oxx[dt] = {}
            #self.oyy[dt] = {}
            #self.txy[dt] = {}
            #self.angle[dt] = {}
            #self.majorP[dt] = {}
            #self.minorP[dt] = {}
            #self.ovmShear[dt] = {}

        #etype = data[0][0]
        #for line in data:
            #if etype == 'CTRIA3':
                #(etype, eid,
                 #fiber_dist1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                 #fiber_dist2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                #self.eType[eid] = etype
                #cen = 0  # CEN/3
                #self.fiberCurvature[eid] = {cen : [fiber_dist1, fiber_dist2]}
                #self.oxx[dt][eid] = {cen : [ox1, ox2]}
                #self.oyy[dt][eid] = {cen : [oy1, oy2]}
                #self.txy[dt][eid] = {cen : [txy1, txy2]}
                #self.angle[dt][eid] = {cen : [angle1, angle2]}
                #self.majorP[dt][eid] = {cen : [o11, o12]}
                #self.minorP[dt][eid] = {cen : [o21, o22]}
                #self.ovmShear[dt][eid] = {cen : [ovm1, ovm2]}
            #elif etype == 'CQUAD4':
                #if len(line) == 18:
                    ## Centroid
                    #(etype, eid,
                     #fiber_dist1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                     #fiber_dist2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                    #nid = 0 # CEN/4
                    #self.eType[eid] = etype
                    #self.fiberCurvature[eid] = {nid : [fiber_dist1, fiber_dist2]}
                    #self.oxx[dt][eid] = {nid : [ox1, ox2]}
                    #self.oyy[dt][eid] = {nid : [oy1, oy2]}
                    #self.txy[dt][eid] = {nid : [txy1, txy2]}
                    #self.angle[dt][eid] = {nid : [angle1, angle2]}
                    #self.majorP[dt][eid] = {nid : [o11, o12]}
                    #self.minorP[dt][eid] = {nid : [o21, o22]}
                    #self.ovmShear[dt][eid] = {nid : [ovm1, ovm2]}
                #elif len(line) == 19:
                    ## Centroid - bilinear
                    #(etype, eid, nid,
                     #fiber_dist1, ox1, oy1, txy1, angle1, o11, o21, ovm1,
                     #fiber_dist2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                    #if nid == 'CEN/4':
                        #nid = 0
                    ##assert isinstance(nid, int), nid
                    #self.eType[eid] = etype
                    #self.fiberCurvature[eid] = {nid: [fiber_dist1, fiber_dist2]}
                    #self.oxx[dt][eid] = {nid: [ox1, ox2]}
                    #self.oyy[dt][eid] = {nid: [oy1, oy2]}
                    #self.txy[dt][eid] = {nid: [txy1, txy2]}
                    #self.angle[dt][eid] = {nid: [angle1, angle2]}
                    #self.majorP[dt][eid] = {nid: [o11, o12]}
                    #self.minorP[dt][eid] = {nid: [o21, o22]}
                    #self.ovmShear[dt][eid] = {nid: [ovm1, ovm2]}
                #elif len(line) == 17:
                    ## Bilinear node
                    #(nid,
                     #fiber_dist1, ox1, oy1, oxy1, angle1, o11, o21, ovm1,
                     #fiber_dist2, ox2, oy2, txy2, angle2, o12, o22, ovm2) = line
                    #assert isinstance(nid, int), nid
                    #self.fiberCurvature[eid][nid] = [fiber_dist1, fiber_dist2]
                    #self.oxx[dt][eid][nid] = [ox1, ox2]
                    #self.oyy[dt][eid][nid] = [oy1, oy2]
                    #self.txy[dt][eid][nid] = [oxy1, txy2]
                    #self.angle[dt][eid][nid] = [angle1, angle2]
                    #self.majorP[dt][eid][nid] = [o11, o12]
                    #self.minorP[dt][eid][nid] = [o21, o22]
                    #self.ovmShear[dt][eid][nid] = [ovm1, ovm2]
                #else:
                    #msg = 'line=%r not supported...len=%i' % (line, len(line))
                    #raise NotImplementedError(msg)
            #else:
                #msg = 'eType=%r not supported...' % etype
                #raise NotImplementedError(msg)

    #def delete_transient(self, dt):
        ##del self.fiberCurvature[dt]
        #del self.oxx[dt]
        #del self.oyy[dt]
        #del self.txy[dt]
        #del self.angle[dt]
        #del self.majorP[dt]
        #del self.minorP[dt]
        #del self.ovmShear[dt]

    #def get_transients(self):
        #k = self.oxx.keys()
        #k.sort()
        #return k

    #def add_new_transient(self, dt):
        #"""
        #initializes the transient variables
        #"""
        #self.dt = dt
        ##self.fiberCurvature[dt] = {}
        #self.oxx[dt] = {}
        #self.oyy[dt] = {}
        #self.txy[dt] = {}
        #self.angle[dt] = {}
        #self.majorP[dt] = {}
        #self.minorP[dt] = {}
        #self.ovmShear[dt] = {}

    #def add_new_eid(self, etype, dt, eid, node_id,
                    #fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        ##if eid in self.oxx:
            ##return self.add(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

        ##assert eid not in self.oxx
        #self.eType[eid] = etype
        #self.fiberCurvature[eid] = {node_id: [fiber_dist]}
        #assert isinstance(eid, int), eid
        #assert isinstance(node_id, int), node_id
        #self.oxx[eid] = {node_id: [oxx]}
        #self.oyy[eid] = {node_id: [oyy]}
        #self.txy[eid] = {node_id: [txy]}
        #self.angle[eid] = {node_id: [angle]}
        #self.majorP[eid] = {node_id: [majorP]}
        #self.minorP[eid] = {node_id: [minorP]}
        #self.ovmShear[eid] = {node_id: [ovm]}
        #msg = "eid=%s node_id=%s fiber_dist=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" % (
            #eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)
        ##if node_id == 0:
            ##raise RuntimeError(msg)

    #def add_new_eid_sort1(self, etype, dt, eid, node_id, fiber_dist,
                          #oxx, oyy, txy, angle, majorP, minorP, ovm):
        ##msg = "dt=%s eid=%s node_id=%s fiber_dist=%g oxx=%g major=%g vm=%g" % (dt, eid, node_id, fiber_dist, oxx, majorP, ovm)
        ##if eid in self.ovmShear[dt]:
        ##    return self.add(eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

        #assert isinstance(eid, int), eid
        #assert isinstance(node_id, int), node_id
        #self.eType[eid] = etype
        #if dt not in self.oxx:
            #self.add_new_transient(dt)
        #self.fiberCurvature[eid] = {node_id: [fiber_dist]}
        #self.oxx[dt][eid] = {node_id: [oxx]}
        #self.oyy[dt][eid] = {node_id: [oyy]}
        #self.txy[dt][eid] = {node_id: [txy]}
        #self.angle[dt][eid] = {node_id: [angle]}
        #self.majorP[dt][eid] = {node_id: [majorP]}
        #self.minorP[dt][eid] = {node_id: [minorP]}
        #self.ovmShear[dt][eid] = {node_id: [ovm]}
        ##if node_id == 0:
            ##msg = "dt=%s eid=%s node_id=%s " % (dt, eid, node_id) #fiber_dist=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(dt,eid,node_id,fiber_dist,oxx,oyy,txy,angle,majorP,minorP,ovm)
            ##raise ValueError(msg)

    #def add(self, dt, eid, node_id,
            #fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        ##msg = "eid=%s node_id=%s fiber_dist=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (
            ##eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #assert isinstance(eid, int), eid
        #assert isinstance(node_id, int), node_id
        #self.fiberCurvature[eid][node_id].append(fiber_dist)
        #self.oxx[eid][node_id].append(oxx)
        #self.oyy[eid][node_id].append(oyy)
        #self.txy[eid][node_id].append(txy)
        #self.angle[eid][node_id].append(angle)
        #self.majorP[eid][node_id].append(majorP)
        #self.minorP[eid][node_id].append(minorP)
        #self.ovmShear[eid][node_id].append(ovm)
        ##if node_id == 0:
            ##raise ValueError(msg)

    #def add_sort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        ##msg = "dt=%s eid=%s node_id=%s fiber_dist=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" % (dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)
        #assert eid is not None, eid
        #assert isinstance(node_id, int), node_id
        #self.fiberCurvature[eid][node_id].append(fiber_dist)
        #self.oxx[dt][eid][node_id].append(oxx)
        #self.oyy[dt][eid][node_id].append(oyy)
        #self.txy[dt][eid][node_id].append(txy)
        #self.angle[dt][eid][node_id].append(angle)
        #self.majorP[dt][eid][node_id].append(majorP)
        #self.minorP[dt][eid][node_id].append(minorP)
        #self.ovmShear[dt][eid][node_id].append(ovm)
        ##if node_id == 0:
            ##raise ValueError(msg)

    #def addNewNode(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert eid is not None, eid
        #assert isinstance(node_id, int), node_id
        ##assert node_id not in self.oxx[eid]
        #self.fiberCurvature[eid][node_id] = [fiber_dist]
        #self.oxx[eid][node_id] = [oxx]
        #self.oyy[eid][node_id] = [oyy]
        #self.txy[eid][node_id] = [txy]
        #self.angle[eid][node_id] = [angle]
        #self.majorP[eid][node_id] = [majorP]
        #self.minorP[eid][node_id] = [minorP]
        #self.ovmShear[eid][node_id] = [ovm]
        ##msg = "eid=%s node_id=%s fiber_dist=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)
        ##print msg
        ##if node_id == 0:
            ##raise ValueError(msg)

    #def addNewNodeSort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert eid is not None
        #assert isinstance(node_id, int), node_id
        ##msg = "dt=%s eid=%s node_id=%s fiber_dist=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)
        ##assert node_id not in self.oxx[dt][eid]
        #self.fiberCurvature[eid][node_id] = [fiber_dist]
        #self.oxx[dt][eid][node_id] = [oxx]
        #self.oyy[dt][eid][node_id] = [oyy]
        #self.txy[dt][eid][node_id] = [txy]
        #self.angle[dt][eid][node_id] = [angle]
        #self.majorP[dt][eid][node_id] = [majorP]
        #self.minorP[dt][eid][node_id] = [minorP]
        #self.ovmShear[dt][eid][node_id] = [ovm]
        ##if node_id == 0:
            ##raise ValueError(msg)

    #def getHeaders(self):
        #if self.is_fiber_distance():
            #headers = ['fiberDist']
        #else:
            #headers = ['curvature']
        #headers += ['oxx', 'oyy', 'txy', 'majorP', 'minorP']
        #if self.is_von_mises():
            #headers.append('oVonMises')
        #else:
            #headers.append('maxShear')
        #return headers

    #def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #if self.nonlinear_factor is not None:
            #return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        #msg_pack, nnodes, cen = _get_plate_msg(self)

        #msg = []
        #msg += header + msg_pack
        #eids = sorted(self.oxx.keys())
        #if self.element_type in [64, 70, 75, 82, 144]:   # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4 bilinear
            #for eid in eids:
                #out = self._write_f06_quad4_bilinear(eid, nnodes, cen)
                #msg.append(out)
        #elif self.element_type in [33, 74]:  # CQUAD4 linear, CTRIA3
            #for eid in eids:
                #out = self._write_f06_tri3(eid)
        #else:
            #raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))

        #msg.append(page_stamp % page_num)
        #f.write(''.join(msg))
        #msg = ['']
        #return page_num

    #def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #msg_pack, nnodes, cen = _get_plate_msg(self)

        #dts = list(self.oxx.keys())
        #dts.sort()
        #if isinstance(dts[0], int):
            #dt_msg = ' %s = %%-10i\n' % self.data_code['name']
        #else:
            #dt_msg = ' %s = %%10.4E\n' % self.data_code['name']

        #msg = []
        #itime = 0
        #ntimes = len(dts)
        #if self.element_type in [64, 70, 75, 82, 144]:  # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4 bilinear
            #for dt, oxx in sorted(iteritems(self.oxx)):
                ##header[1] = dt_msg % dt
                #header = _eigenvalue_header(self, header, itime, ntimes, dt)
                #msg.append('\n'.join(header + msg_pack))
                #for eid in sorted(oxx.keys()):
                    #out = self._write_f06_quad4_bilinear_transient(dt, eid, nnodes, cen)
                    #msg.append(out)
                #msg.append(page_stamp % page_num)
                #f.write(''.join(msg))
                #msg = ['']
                #page_num += 1
                #itime += 1
        #elif self.element_type in [33, 74]:
            #for dt, oxx in sorted(iteritems(self.oxx)):
                ##header[1] = dt_msg % dt
                #header = _eigenvalue_header(self, header, itime, ntimes, dt)
                #msg.append('\n'.join(header + msg_pack))
                #for eid in sorted(oxx.keys()):
                    #out = self._write_f06_tri3_transient(dt, eid)
                    #msg.append(out)
                #msg.append(page_stamp % page_num)
                #f.write(''.join(msg))
                #msg = ['']
                #page_num += 1
                #itime += 1
        #else:
            #raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))
        #return page_num - 1

    #def _write_f06_quad4_bilinear(self, eid, n, cen):
        #"""Writes the static CQUAD4 bilinear"""
        #msg = ['']
        #nids = sorted(self.oxx[eid].keys())
        #for nid in nids:
            #for ilayer in range(len(self.oxx[eid][nid])):
                #fiber_dist = self.fiberCurvature[eid][nid][ilayer]
                #oxx = self.oxx[eid][nid][ilayer]
                #oyy = self.oyy[eid][nid][ilayer]
                #txy = self.txy[eid][nid][ilayer]
                #angle = self.angle[eid][nid][ilayer]
                #major = self.majorP[eid][nid][ilayer]
                #minor = self.minorP[eid][nid][ilayer]
                #ovm = self.ovmShear[eid][nid][ilayer]
                #([fiber_dist, oxx, oyy, txy, major, minor, ovm], is_all_zeros) = writeFloats13E([fiber_dist, oxx, oyy, txy, major, minor, ovm])
                #([angle], is_all_zeros) = writeFloats8p4F([angle])

                #if nid == 0 and ilayer == 0:
                    #msg.append('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % (eid, cen, fiber_dist, oxx, oyy, txy, angle, major, minor, ovm))
                #elif ilayer == 0:
                    #msg.append('   %8s %8i  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % ('', nid, fiber_dist, oxx, oyy, txy, angle, major, minor, ovm))
                #elif ilayer == 1:
                    #msg.append('   %8s %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n\n' % ('', '', fiber_dist, oxx, oyy, txy, angle, major, minor, ovm))
                #else:
                    #raise Exception('Invalid option for cquad4')
        #return ''.join(msg)

    #def _write_f06_quad4_bilinear_transient(self, dt, eid, n, cen):
        #"""Writes the transient CQUAD4 bilinear"""
        #msg = ['']
        #nids = sorted(self.oxx[dt][eid].keys())
        #for nid in nids:
            #for ilayer in range(len(self.oxx[dt][eid][nid])):
                #fiber_dist = self.fiberCurvature[eid][nid][ilayer]
                #oxx = self.oxx[dt][eid][nid][ilayer]
                #oyy = self.oyy[dt][eid][nid][ilayer]
                #txy = self.txy[dt][eid][nid][ilayer]
                #angle = self.angle[dt][eid][nid][ilayer]
                #major = self.majorP[dt][eid][nid][ilayer]
                #minor = self.minorP[dt][eid][nid][ilayer]
                #ovm = self.ovmShear[dt][eid][nid][ilayer]
                #([fiber_dist, oxx, oyy, txy, major, minor, ovm], is_all_zeros) = writeFloats13E([fiber_dist, oxx, oyy, txy, major, minor, ovm])
                #([angle], is_all_zeros) = writeFloats8p4F([angle])

                #if nid == 0 and ilayer == 0:
                    #msg.append('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % (eid, cen, fiber_dist, oxx, oyy, txy, angle, major, minor, ovm))
                #elif ilayer == 0:
                    #msg.append('   %8s %8i  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % ('', nid, fiber_dist, oxx, oyy, txy, angle, major, minor, ovm))
                #elif ilayer == 1:
                    #msg.append('   %8s %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n\n' % ('', '', fiber_dist, oxx, oyy, txy, angle, major, minor, ovm))
                #else:
                    #raise RuntimeError('Invalid option for cquad4')
        #return ''.join(msg)

    #def _write_f06_tri3(self, eid):
        #msg = ['']
        #oxx_nodes = self.oxx[eid].keys()
        #for nid in sorted(oxx_nodes):
            #for ilayer in range(len(self.oxx[eid][nid])):
                #fiber_dist = self.fiberCurvature[eid][nid][ilayer]
                #oxx = self.oxx[eid][nid][ilayer]
                #oyy = self.oyy[eid][nid][ilayer]
                #txy = self.txy[eid][nid][ilayer]
                #angle = self.angle[eid][nid][ilayer]
                #major = self.majorP[eid][nid][ilayer]
                #minor = self.minorP[eid][nid][ilayer]
                #ovm = self.ovmShear[eid][nid][ilayer]
                #([fiber_dist, oxx, oyy, txy, major, minor, ovm], is_all_zeros) = writeFloats13E([fiber_dist, oxx, oyy, txy, major, minor, ovm])
                #([angle], is_all_zeros) = writeFloats8p4F([angle])

                #if ilayer == 0:
                    #msg.append('0  %6i   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %-s\n' % (eid, fiber_dist, oxx, oyy, txy, angle, major, minor, ovm.rstrip()))
                #else:
                    #msg.append('   %6s   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %-s\n' % ('', fiber_dist, oxx, oyy, txy, angle, major, minor, ovm.rstrip()))
        #return ''.join(msg)

    #def _write_f06_tri3_transient(self, dt, eid):
        #msg = ['']
        #oxx_nodes = self.oxx[dt][eid].keys()
        #for nid in sorted(oxx_nodes):
            #for ilayer in range(len(self.oxx[dt][eid][nid])):
                #fiber_dist = self.fiberCurvature[eid][nid][ilayer]
                #oxx = self.oxx[dt][eid][nid][ilayer]
                #oyy = self.oyy[dt][eid][nid][ilayer]
                #txy = self.txy[dt][eid][nid][ilayer]
                #angle = self.angle[dt][eid][nid][ilayer]
                #major = self.majorP[dt][eid][nid][ilayer]
                #minor = self.minorP[dt][eid][nid][ilayer]
                #ovm = self.ovmShear[dt][eid][nid][ilayer]
                #([fiber_dist, oxx, oyy, txy, major, minor, ovm], is_all_zeros) = writeFloats13E([fiber_dist, oxx, oyy, txy, major, minor, ovm])
                #([angle], is_all_zeros) = writeFloats8p4F([angle])

                #if ilayer == 0:
                    #msg.append('0  %6i   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %s\n' % (eid, fiber_dist, oxx, oyy, txy, angle, major, minor, ovm))
                #else:
                    #msg.append('   %6s   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %s\n' % ('', fiber_dist, oxx, oyy, txy, angle, major, minor, ovm))
        #return ''.join(msg)


#class RealPlateStrain(StrainObject):
    #"""
    #::

      ## ??? - is this just 11
      #ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)
        #ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES

      ## s_code=11
                             #S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      #ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)
        #ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       VON MISES

      ## s_code=15
                             #S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
      #ELEMENT      FIBER                STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)
        #ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES

      ## s_code=10
                             #S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN
      #ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)          MAX
        #ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR
    #"""
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StrainObject.__init__(self, data_code, isubcase)
        #self.eType = {}

        #self.code = [self.format_code, self.sort_code, self.s_code]

        ##print self.data_code
        #self.fiberCurvature = {}
        #self.exx = {}
        #self.eyy = {}
        #self.exy = {}
        #self.angle = {}
        #self.majorP = {}
        #self.minorP = {}
        #self.evmShear = {}

        #self.dt = dt
        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
                #self.addNewNode = self.addNewNodeSort1
        #else:
            #raise NotImplementedError('SORT2')

    #def get_stats(self):
        #nelements = len(self.eType)

        #msg = self.get_data_code()
        #if self.nonlinear_factor is not None:  # transient
            #ntimes = len(self.exx)
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
        #msg.append('  eType, fiberCurvature, exx, eyy, exy, angle, '
                   #'majorP, minorP, evmShear\n')
        #msg.append('  scode=%s stress_bits=%s is_stress=%s dist=%s vm=%s\n' % (
            #self.s_code, self.stress_bits, self.is_stress(), self.is_fiber_distance(), self.is_von_mises()))
        #return msg

    #def add_f06_data(self, data, transient, data_code=None):
        #if transient is None:
            #etype = data[0][0]
            #for line in data:
                #if etype == 'CTRIA3':
                    #cen = 0
                    #(etype, eid,
                     #fiber_dist1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                     #fiber_dist2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    #self.eType[eid] = etype
                    #self.fiberCurvature[eid] = {cen : [fiber_dist1, fiber_dist2]}
                    #self.exx[eid] = {cen : [ex1, ex2]}
                    #self.eyy[eid] = {cen : [ey1, ey2]}
                    #self.exy[eid] = {cen : [exy1, exy2]}
                    #self.angle[eid] = {cen : [angle1, angle2]}
                    #self.majorP[eid] = {cen : [e11, e12]}
                    #self.minorP[eid] = {cen : [e21, e22]}
                    #self.evmShear[eid] = {cen : [evm1, evm2]}
                #elif etype == 'CQUAD4':
                    #if len(line) == 19:  # Centroid - bilinear
                        #(etype, eid, nid,
                         #fiber_dist1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                         #fiber_dist2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        #if nid == 'CEN/4':
                            #nid = 0
                        ##assert isinstance(nid, int), nid
                        #self.eType[eid] = etype
                        #self.fiberCurvature[eid] = {nid : [fiber_dist1, fiber_dist2]}
                        #self.exx[eid] = {nid : [ex1, ex2]}
                        #self.eyy[eid] = {nid : [ey1, ey2]}
                        #self.exy[eid] = {nid : [exy1, exy2]}
                        #self.angle[eid] = {nid : [angle1, angle2]}
                        #self.majorP[eid] = {nid : [e11, e12]}
                        #self.minorP[eid] = {nid : [e21, e22]}
                        #self.evmShear[eid] = {nid : [evm1, evm2]}
                    #elif len(line) == 18:  # Centroid
                        #(etype, eid,
                         #fiber_dist1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                         #fiber_dist2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        #nid = 0 # CEN/4
                        ##assert isinstance(nid, cen), nid
                        #self.eType[eid] = etype
                        #self.fiberCurvature[eid] = {nid : [fiber_dist1, fiber_dist2]}
                        #self.exx[eid] = {nid : [ex1, ex2]}
                        #self.eyy[eid] = {nid : [ey1, ey2]}
                        #self.exy[eid] = {nid : [exy1, exy2]}
                        #self.angle[eid] = {nid : [angle1, angle2]}
                        #self.majorP[eid] = {nid : [e11, e12]}
                        #self.minorP[eid] = {nid : [e21, e22]}
                        #self.evmShear[eid] = {nid : [evm1, evm2]}
                    #elif len(line) == 17:  # Bilinear node
                        #(nid,
                         #fiber_dist1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                         #fiber_dist2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                        #assert isinstance(nid, int), nid
                        #self.fiberCurvature[eid][nid] = [fiber_dist1, fiber_dist2]
                        #self.exx[eid][nid] = [ex1, ex2]
                        #self.eyy[eid][nid] = [ey1, ey2]
                        #self.exy[eid][nid] = [exy1, exy2]
                        #self.angle[eid][nid] = [angle1, angle2]
                        #self.majorP[eid][nid] = [e11, e12]
                        #self.minorP[eid][nid] = [e21, e22]
                        #self.evmShear[eid][nid] = [evm1, evm2]
                    #else:
                        ##assert len(line) == 19, 'len(line)=%s' % len(line)
                        #msg = 'line=%r not supported...len=%i' % (line, len(line))
                        #raise NotImplementedError(msg)
                #else:
                    #msg = 'line=%s not supported...' % line
                    #raise NotImplementedError(msg)
            #return
        #etype = data[0][0]
        #assert 'name' in self.data_code, self.data_code

        #dt = transient[1]
        #if dt not in self.exx:
            #self.exx[dt] = {}
            #self.eyy[dt] = {}
            #self.exy[dt] = {}
            #self.angle[dt] = {}
            #self.majorP[dt] = {}
            #self.minorP[dt] = {}
            #self.evmShear[dt] = {}

        #for line in data:
            #etype = data[0][0]
            #if etype == 'CTRIA3':
                #cen = 0 #  CEN/3
                #(etype, eid,
                 #fiber_dist1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                 #fiber_dist2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                #self.eType[eid] = etype
                #self.fiberCurvature[eid] = {cen : [fiber_dist1, fiber_dist2]}
                #self.exx[dt][eid] = {cen : [ex1, ex2]}
                #self.eyy[dt][eid] = {cen : [ey1, ey2]}
                #self.exy[dt][eid] = {cen : [exy1, exy2]}
                #self.angle[dt][eid] = {cen : [angle1, angle2]}
                #self.majorP[dt][eid] = {cen : [e11, e12]}
                #self.minorP[dt][eid] = {cen : [e21, e22]}
                #self.evmShear[dt][eid] = {cen : [evm1, evm2]}
            #elif etype == 'CQUAD4':
                #if len(line) == 19:  # Centroid - bilinear
                    #(etype, eid, nid,
                     #fiber_dist1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                     #fiber_dist2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    ##assert isinstance(nid, int), nid
                    #if nid == 'CEN/4':
                        #nid = 0

                    #self.eType[eid] = etype
                    #self.fiberCurvature[eid] = {nid: [fiber_dist1, fiber_dist2]}
                    #self.exx[dt][eid] = {nid: [ex1, ex2]}
                    #self.eyy[dt][eid] = {nid: [ey1, ey2]}
                    #self.exy[dt][eid] = {nid: [exy1, exy2]}
                    #self.angle[dt][eid] = {nid: [angle1, angle2]}
                    #self.majorP[dt][eid] = {nid: [e11, e12]}
                    #self.minorP[dt][eid] = {nid: [e21, e22]}
                    #self.evmShear[dt][eid] = {nid: [evm1, evm2]}
                #elif len(line) == 18:  # Centroid
                    #(etype, eid,
                     #fiber_dist1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                     #fiber_dist2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    #self.eType[eid] = etype
                    ##assert isinstance(nid, int), nid
                    #nid = 0
                    #self.fiberCurvature[eid] = {nid : [fiber_dist1, fiber_dist2]}
                    #self.exx[dt][eid] = {nid : [ex1, ex2]}
                    #self.eyy[dt][eid] = {nid : [ey1, ey2]}
                    #self.exy[dt][eid] = {nid : [exy1, exy2]}
                    #self.angle[dt][eid] = {nid : [angle1, angle2]}
                    #self.majorP[dt][eid] = {nid : [e11, e12]}
                    #self.minorP[dt][eid] = {nid : [e21, e22]}
                    #self.evmShear[dt][eid] = {nid : [evm1, evm2]}
                #elif len(line) == 17:  # Bilinear node
                    #(nid,
                     #fiber_dist1, ex1, ey1, exy1, angle1, e11, e21, evm1,
                     #fiber_dist2, ex2, ey2, exy2, angle2, e12, e22, evm2) = line
                    #assert isinstance(nid, int), nid
                    #self.fiberCurvature[eid][nid] = [fiber_dist1, fiber_dist2]
                    #self.exx[dt][eid][nid] = [ex1, ex2]
                    #self.eyy[dt][eid][nid] = [ey1, ey2]
                    #self.exy[dt][eid][nid] = [exy1, exy2]
                    #self.angle[dt][eid][nid] = [angle1, angle2]
                    #self.majorP[dt][eid][nid] = [e11, e12]
                    #self.minorP[dt][eid][nid] = [e21, e22]
                    #self.evmShear[dt][eid][nid] = [evm1, evm2]
                #else:
                    #msg = 'line=%r not supported...len=%i' % (line, len(line))
                    #raise NotImplementedError(msg)
            #else:
                #msg = 'eType=%r is not supported...' % etype
                #raise NotImplementedError(msg)

    #def delete_transient(self, dt):
        ##del self.fiberCurvature[dt]
        #del self.exx[dt]
        #del self.eyy[dt]
        #del self.exy[dt]
        #del self.angle[dt]
        #del self.majorP[dt]
        #del self.minorP[dt]
        #del self.evmShear[dt]

    #def get_transients(self):
        #k = self.exx.keys()
        #k.sort()
        #return k

    #def add_new_transient(self, dt):
        #"""
        #initializes the transient variables
        #"""
        #self.exx[dt] = {}
        #self.eyy[dt] = {}
        #self.exy[dt] = {}
        #self.angle[dt] = {}
        #self.majorP[dt] = {}
        #self.minorP[dt] = {}
        #self.evmShear[dt] = {}

    #def add_new_eid(self, etype, dt, eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        ##msg = "eid=%s node_id=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #assert isinstance(node_id, int), node_id

        ##if node_id != 'C':  # centroid
            ##assert 0 < node_id < 1000000000, 'node_id=%s %s' % (node_id, msg)
        ##assert eid not in self.exx
        #self.eType[eid] = etype
        #assert isinstance(node_id, int), node_id
        #self.fiberCurvature[eid] = {node_id : [curvature]}
        #self.exx[eid] = {node_id : [exx]}
        #self.eyy[eid] = {node_id : [eyy]}
        #self.exy[eid] = {node_id : [exy]}
        #self.angle[eid] = {node_id : [angle]}
        #self.majorP[eid] = {node_id : [majorP]}
        #self.minorP[eid] = {node_id : [minorP]}
        #self.evmShear[eid] = {node_id : [evm]}

    #def add_new_eid_sort1(self, etype, dt, eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        ##msg = "eid=%s node_id=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #assert isinstance(node_id, int), node_id
        ##print msg

        ##if node_id != 'C':  # centroid
            ##assert 0 < node_id < 1000000000, 'node_id=%s %s' % (node_id, msg)

        #if dt not in self.exx:
            #self.add_new_transient(dt)
        ##if eid in self.evmShear[dt]:  # SOL200, erase the old result
            ##nid = node_id
            ##msg = "dt=%s eid=%s node_id=%s fiber_dist=%s oxx=%s major=%s vm=%s" %(dt, eid, node_id, str(self.fiberCurvature[eid][nid]), str(self.oxx[dt][eid][nid]),str(self.majorP[dt][eid][nid]),str(self.ovmShear[dt][eid][nid]))
            ##self.delete_transient(dt)
            ##self.add_new_transient()

        #self.eType[eid] = etype
        #assert isinstance(node_id, int), node_id
        #self.fiberCurvature[eid] = {node_id: [curvature]}
        #self.exx[dt][eid] = {node_id: [exx]}
        #self.eyy[dt][eid] = {node_id: [eyy]}
        #self.exy[dt][eid] = {node_id: [exy]}
        #self.angle[dt][eid] = {node_id: [angle]}
        #self.majorP[dt][eid] = {node_id: [majorP]}
        #self.minorP[dt][eid] = {node_id: [minorP]}
        #self.evmShear[dt][eid] = {node_id: [evm]}

    #def add(self, dt, eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        ##msg = "eid=%s node_id=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        ##if node_id != 'C':  # centroid
            ##assert 0 < node_id < 1000000000, 'node_id=%s' % node_id
        #assert isinstance(node_id, int), node_id
        #self.fiberCurvature[eid][node_id].append(curvature)
        #self.exx[eid][node_id].append(exx)
        #self.eyy[eid][node_id].append(eyy)
        #self.exy[eid][node_id].append(exy)
        #self.angle[eid][node_id].append(angle)
        #self.majorP[eid][node_id].append(majorP)
        #self.minorP[eid][node_id].append(minorP)
        #self.evmShear[eid][node_id].append(evm)

    #def add_sort1(self, dt, eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        ##msg = "eid=%s node_id=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        ##if node_id != 'C':  # centroid
            ##assert 0 < node_id < 1000000000, 'node_id=%s' % (node_id)

        #assert isinstance(node_id, int), node_id
        #self.fiberCurvature[eid][node_id].append(curvature)
        #self.exx[dt][eid][node_id].append(exx)
        #self.eyy[dt][eid][node_id].append(eyy)
        #self.exy[dt][eid][node_id].append(exy)
        #self.angle[dt][eid][node_id].append(angle)
        #self.majorP[dt][eid][node_id].append(majorP)
        #self.minorP[dt][eid][node_id].append(minorP)
        #self.evmShear[dt][eid][node_id].append(evm)

    #def addNewNode(self, dt, eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #msg = "eid=%s node_id=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        #assert node_id not in self.exx[eid], msg
        #assert isinstance(node_id, int), node_id

        #self.fiberCurvature[eid][node_id] = [curvature]
        #self.exx[eid][node_id] = [exx]
        #self.eyy[eid][node_id] = [eyy]
        #self.exy[eid][node_id] = [exy]
        #self.angle[eid][node_id] = [angle]
        #self.majorP[eid][node_id] = [majorP]
        #self.minorP[eid][node_id] = [minorP]
        #self.evmShear[eid][node_id] = [evm]

    #def addNewNodeSort1(self, dt, eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm):
        #msg = "eid=%s node_id=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" % (eid, node_id, curvature, exx, eyy, exy, angle, majorP, minorP, evm)
        ##assert node_id not in self.exx[eid], msg
        #assert isinstance(node_id, int), node_id
        #self.fiberCurvature[eid][node_id] = [curvature]
        #self.exx[dt][eid][node_id] = [exx]
        #self.eyy[dt][eid][node_id] = [eyy]
        #self.exy[dt][eid][node_id] = [exy]
        #self.angle[dt][eid][node_id] = [angle]
        #self.majorP[dt][eid][node_id] = [majorP]
        #self.minorP[dt][eid][node_id] = [minorP]
        #self.evmShear[dt][eid][node_id] = [evm]

    #def getHeaders(self):
        #if self.is_fiber_distance():
            #headers = ['fiberDist']
        #else:
            #headers = ['curvature']

        #headers += ['exx', 'eyy', 'exy', 'eMajor', 'eMinor']
        #if self.is_von_mises():
            #headers.append('eVonMises')
        #else:
            #headers.append('maxShear')
        #return headers

    #def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #if self.nonlinear_factor is not None:
            #return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        #msg_pack, nnodes, cen = _get_plate_msg(self)

        #msg = []
        #msg += header + msg_pack
        #eids = sorted(self.exx.keys())

        #if self.element_type in [64, 70, 75, 82, 144]:
            #for eid in eids:
                #out = self._write_f06_quad4_bilinear(eid, nnodes, cen)
                #msg.append(out)
        #elif self.element_type in [33, 74]:
            #for eid in eids:
                #out = self._write_f06_tri3(eid)
                #msg.append(out)
        #else:
            #raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))

        #msg.append(page_stamp % page_num)
        #f.write(''.join(msg))
        #return page_num

    #def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #msg_pack, nnodes, cen = _get_plate_msg(self)

        #dts = list(self.exx.keys())
        #dts.sort()
        #if isinstance(dts[0], int):
            #dt_msg = ' %s = %%-10i\n' % self.data_code['name']
        #else:
            #dt_msg = ' %s = %%10.4E\n' % self.data_code['name']

        #msg = []
        #itime = 0
        #ntimes = len(dts)
        #if self.element_type in [64, 70, 75, 82, 144]:
            #for dt, oxx in sorted(iteritems(self.exx)):
                #header = _eigenvalue_header(self, header, itime, ntimes, dt)
                ##header[1] = dt_msg % dt
                #msg.append('\n'.join(header + msg_pack))
                #for eid in sorted(oxx.keys()):
                    #out = self._write_f06_quad4_bilinear_transient(dt, eid, nnodes, cen)
                    #msg.append(out)
                #msg.append(page_stamp % page_num)
                #page_num += 1
                #itime += 1
        #elif self.element_type in [33, 74]:
            #for dt, oxx in sorted(iteritems(self.exx)):
                #header = _eigenvalue_header(self, header, itime, ntimes, dt)
                ##header[1] = dt_msg % dt
                #msg.append('\n'.join(header + msg_pack))
                #for eid in sorted(oxx.keys()):
                    #out = self._write_f06_tri3_transient(dt, eid)
                    #msg.append(out)
                #msg.append(page_stamp % page_num)
                #page_num += 1
                #itime += 1
        #else:
            #raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))
        #f.write(''.join(msg))
        #return page_num - 1

    #def _write_f06_quad4_bilinear(self, eid, n, cen):
        #"""Writes the static CQUAD4 bilinear"""
        #msg = ['']
        #nids = sorted(self.exx[eid].keys())
        #for nid in nids:
            #for ilayer in range(len(self.exx[eid][nid])):
                #fiber_dist = self.fiberCurvature[eid][nid][ilayer]
                #exx = self.exx[eid][nid][ilayer]
                #eyy = self.eyy[eid][nid][ilayer]
                #exy = self.exy[eid][nid][ilayer]
                #angle = self.angle[eid][nid][ilayer]
                #major = self.majorP[eid][nid][ilayer]
                #minor = self.minorP[eid][nid][ilayer]
                #evm = self.evmShear[eid][nid][ilayer]
                #([fiber_dist, exx, eyy, exy, major, minor, evm], is_all_zeros) = writeFloats13E([fiber_dist, exx, eyy, exy, major, minor, evm])
                #([angle], is_all_zeros) = writeFloats8p4F([angle])

                #if nid == 0 and ilayer == 0:
                    #msg.append('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % (eid, cen, fiber_dist, exx, eyy, exy, angle, major, minor, evm))
                #elif ilayer == 0:
                    #msg.append('   %8s %8i  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n' % ('', nid, fiber_dist, exx, eyy, exy, angle, major, minor, evm))
                #elif ilayer == 1:
                    #msg.append('   %8s %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %s\n\n' % ('', '', fiber_dist, exx, eyy, exy, angle, major, minor, evm))
                #else:
                    #raise RuntimeError('Invalid option for cquad4')
        #return ''.join(msg)

    #def _write_f06_quad4_bilinear_transient(self, dt, eid, n, cen):
        #"""Writes the transient CQUAD4 bilinear"""
        #msg = ['']
        #nids = sorted(self.exx[dt][eid].keys())
        #for nid in nids:
            #for ilayer in range(len(self.exx[dt][eid][nid])):
                #fiber_dist = self.fiberCurvature[eid][nid][ilayer]
                #exx = self.exx[dt][eid][nid][ilayer]
                #eyy = self.eyy[dt][eid][nid][ilayer]
                #exy = self.exy[dt][eid][nid][ilayer]
                #angle = self.angle[dt][eid][nid][ilayer]
                #major = self.majorP[dt][eid][nid][ilayer]
                #minor = self.minorP[dt][eid][nid][ilayer]
                #evm = self.evmShear[dt][eid][nid][ilayer]

                #([fiber_dist, exx, eyy, exy, major, minor, evm], is_all_zeros) = writeFloats13E([fiber_dist, exx, eyy, exy, major, minor, evm])
                #([angle], is_all_zeros) = writeFloats8p4F([angle])

                #if nid == 0 and ilayer == 0:
                    #msg.append('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %-s\n' % (eid, cen, fiber_dist, exx, eyy, exy, angle, major, minor, evm.rstrip()))
                #elif ilayer == 0:
                    #msg.append('   %8s %8i  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %-s\n' % ('', nid, fiber_dist, exx, eyy, exy, angle, major, minor, evm.rstrip()))
                #elif ilayer == 1:
                    #msg.append('   %8s %8s  %-13s  %-13s %-13s %-13s   %8s  %-13s %-13s %-s\n\n' % ('', '', fiber_dist, exx, eyy, exy, angle, major, minor, evm.rstrip()))
                #else:
                    #raise RuntimeError('Invalid option for cquad4')
        #return ''.join(msg)

    #def _write_f06_tri3(self, eid):
        #"""Writes the static CTRIA3/CQUAD4 linear"""
        #msg = ['']
        #for nid in sorted(self.exx[eid].keys()):
            #for ilayer in range(len(self.exx[eid][nid])):
                #fiber_dist = self.fiberCurvature[eid][nid][ilayer]
                #exx = self.exx[eid][nid][ilayer]
                #eyy = self.eyy[eid][nid][ilayer]
                #exy = self.exy[eid][nid][ilayer]
                #angle = self.angle[eid][nid][ilayer]
                #major = self.majorP[eid][nid][ilayer]
                #minor = self.minorP[eid][nid][ilayer]
                #evm = self.evmShear[eid][nid][ilayer]

                #([fiber_dist, exx, eyy, exy, major, minor, evm], is_all_zeros) = writeFloats13E([fiber_dist, exx, eyy, exy, major, minor, evm])
                #([angle], is_all_zeros) = writeFloats8p4F([angle])
                #if ilayer == 0:
                    #msg.append('0  %6i   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %-s\n' % (eid, fiber_dist, exx, eyy, exy, angle, major, minor, evm.rstrip()))
                #else:
                    #msg.append('   %6s   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %-s\n' % ('', fiber_dist, exx, eyy, exy, angle, major, minor, evm.rstrip()))
        #return ''.join(msg)

    #def _write_f06_tri3_transient(self, dt, eid):
        #"""Writes the transient CTRIA3/CQUAD4 linear"""
        #msg = ['']
        #for nid in sorted(self.exx[dt][eid]):
            #for ilayer in range(len(self.exx[dt][eid][nid])):
                #fiber_dist = self.fiberCurvature[eid][nid][ilayer]
                #exx = self.exx[dt][eid][nid][ilayer]
                #eyy = self.eyy[dt][eid][nid][ilayer]
                #exy = self.exy[dt][eid][nid][ilayer]
                #angle = self.angle[dt][eid][nid][ilayer]
                #major = self.majorP[dt][eid][nid][ilayer]
                #minor = self.minorP[dt][eid][nid][ilayer]
                #evm = self.evmShear[dt][eid][nid][ilayer]

                #([fiber_dist, exx, eyy, exy, major, minor, evm], is_all_zeros) = writeFloats13E([fiber_dist, exx, eyy, exy, major, minor, evm])
                #([angle], is_all_zeros) = writeFloats8p4F([angle])
                #if ilayer == 0:
                    #msg.append(
                        #'0  %6i   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %-s\n' % (
                            #eid, fiber_dist, exx, eyy, exy, angle, major, minor, evm))
                #else:
                    #msg.append('   %6s   %-13s     %-13s  %-13s  %-13s   %8s   %-13s   %-13s  %-s\n' % (
                        #'', fiber_dist, exx, eyy, exy, angle, major, minor, evm))
        #return ''.join(msg)
