# coding: utf-8
#pylint disable=C0103
from typing import List
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object


class RealCPLSTRNPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        if not is_sort1:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def is_bilinear(self):
        #if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
            #return False
        #elif self.element_type in [144, 64, 82, 70, 75]:  # CQUAD4
            #return True
        raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))

    def build(self):
        """sizes the vectorized attributes of the RealCPLSTRNPlateArray"""
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 275:  # CPLSTS3
            nnodes_per_element = 1
        elif self.element_type == 276:  # CPLSTS4
            nnodes_per_element = 5
        #elif self.element_type == 64:  # CQUAD8
            #nnodes_per_element = 5
        #elif self.element_type == 82:  # CQUADR
            #nnodes_per_element = 5
        #elif self.element_type == 70:  # CTRIAR
            #nnodes_per_element = 4
        #elif self.element_type == 75:  # CTRIA6
            #nnodes_per_element = 4
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
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.element = np.zeros(self.ntotal, dtype='int32')

        #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        self.data = np.zeros((self.ntimes, self.ntotal, 5), dtype='float32')

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fiber_dist1, oxx1, oyy1, txy1, angle1, major_p1, minor_p1, ovm1) = t1
                    (fiber_dist2, oxx2, oyy2, txy2, angle2, major_p2, minor_p2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            fiber_dist1, oxx1, oyy1, txy1, angle1, major_p1, minor_p1, ovm1,
                            fiber_dist2, oxx2, oyy2, txy2, angle2, major_p2, minor_p2, ovm2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def get_stats(self, short=False) -> List[str]:
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
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msgi = '  type=%s ntimes=%i nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, nlayers, ntotal)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, nelements, nnodes, nlayers, ntotal)
            ntimes_word = '1'
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 str(', '.join(headers))))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  data.shape=%s\n' % str(self.data.shape))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    #def get_element_index(self, eids):
        ## elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = np.ravel([np.searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ##ind = np.searchsorted(eids, self.element)
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    #def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  #page_num=1, is_mag_phase=False, is_sort1=True):
        #if header is None:
            #header = []
        #msg, nnodes, cen = _get_plate_msg(self)

        ## write the f06
        #ntimes = self.data.shape[0]

        #eids = self.element_node[:, 0]
        #nids = self.element_node[:, 1]

        ##cen_word = 'CEN/%i' % nnodes
        #cen_word = cen
        #for itime in range(ntimes):
            #dt = self._times[itime]
            #header = _eigenvalue_header(self, header, itime, ntimes, dt)
            #f06_file.write(''.join(header + msg))

            ##print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            ##[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
            #fiber_dist = self.data[itime, :, 0]
            #oxx = self.data[itime, :, 1]
            #oyy = self.data[itime, :, 2]
            #txy = self.data[itime, :, 3]
            #angle = self.data[itime, :, 4]
            #majorP = self.data[itime, :, 5]
            #minorP = self.data[itime, :, 6]
            #ovm = self.data[itime, :, 7]

            #for (i, eid, nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi) in zip(
                 #count(), eids, nids, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
                #[fdi, oxxi, oyyi, txyi, major, minor, ovmi] = write_floats_13e(
                #[fdi, oxxi, oyyi, txyi, major, minor, ovmi])
                #ilayer = i % 2
                ## tria3
                #if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
                    #if ilayer == 0:
                        #f06_file.write('0  %6i   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (
                            #eid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    #else:
                        #f06_file.write('   %6s   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (
                            #'', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))

                #elif self.element_type in [64, 70, 75, 82, 144]:  # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4
                    ## bilinear
                    #if nid == 0 and ilayer == 0:  # CEN
                        #f06_file.write('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (
                            #eid, cen_word, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    #elif ilayer == 0:
                        #f06_file.write('   %8s %8i  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (
                            #'', nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    #elif ilayer == 1:
                        #f06_file.write('   %8s %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n\n' % (
                            #'', '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                #else:
                    #raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))

            #f06_file.write(page_stamp % page_num)
            #page_num += 1
        #return page_num - 1

class RealCPLSTRNPlateStressArray(RealCPLSTRNPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCPLSTRNPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['oxx', 'oyy', 'txy', 'von_mises']
        return headers


class RealCPLSTRNPlateStrainArray(RealCPLSTRNPlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCPLSTRNPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['exx', 'eyy', 'exy', 'von_mises']
        return headers


def _get_clpstsn_msg(self):
    if self.is_von_mises:
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    if self.is_stress:
        #if self.is_fiber_distance:
            #quad_msg_temp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                             #'      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            #tri_msg_temp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                            #'    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        #else:
            #quad_msg_temp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                             #'      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
        tri_msg_temp = ['    ELEMENT\n',
                        '      ID           NORMAL-X          NORMAL-Y          NORMAL-Z           SHEAR           VON MISES\n']



        #'0        26     -4.183719E-05      4.101315E-03      3.378940E-04      1.064501E-03      2.715101E-03'

        #is_bilinear = False
        #cpltsts3_msg = ['                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( C P L S T S 3 )\n'] + tri_msg_temp
        #cquad8_msg = ['                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( C P L S T S 3 )\n'] + tri_msg_temp
        #cquadr_msg = ['                S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( C P L S T S 3 )\n'] + tri_msg_temp

        #is_bilinear = True
        #cquadr_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quad_msg_temp
        #cquad4_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quad_msg_temp

        #is_bilinear = False
        ctria3_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + tri_msg_temp
        ctria6_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + tri_msg_temp
        ctriar_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + tri_msg_temp
    else:
        if self.is_fiber_distance:
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
