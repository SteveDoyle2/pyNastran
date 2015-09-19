from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip, range
from numpy import zeros, searchsorted, unique, ravel

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats12E, _eigenvalue_header


class RealCompositePlateArray(OES_Object):
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
                pass
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def build(self):
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        if self.element_type == 95:  # CQUAD4
            nnodes_per_element = 1
        elif self.element_type == 96:  # CQUAD8
            nnodes_per_element = 1
        elif self.element_type == 97:  # CTRIA3
            nnodes_per_element = 1
        elif self.element_type == 98:  # CTRIA6
            nnodes_per_element = 1
        else:
            raise NotImplementedError('element_name=%s element_type=%s' %(self.element_name, self.element_type))

        self.nnodes = nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True

        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)

        self.element_layer = zeros((self.ntotal, 2), dtype='int32')

        #[o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
        self.data = zeros((self.ntimes, self.ntotal, 9), dtype='float32')

    def add_new_eid(self, eType, dt, eid, layer, o11, o22, t12, t1z, t2z,
                    angle, major, minor, ovm):
        self.add_new_eid_sort1(eType, dt, eid, layer, o11, o22, t12, t1z, t2z,
                               angle, major, minor, ovm)

    def add_new_eid_sort1(self, etype, dt, eid, layer, o11, o22, t12, t1z, t2z,
                          angle, major, minor, ovm):
        self._times[self.itime] = dt
        self.element_layer[self.itotal, :] = [eid, layer]
        self.data[self.itime, self.itotal, :] = [o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
        self.itotal += 1
        self.ielement += 1

    def add(self, dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                 major, minor, ovm):
        self.add_sort1(dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                       major, minor, ovm)

    def add_sort1(self, dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                  major, minor, ovm):
        assert eid is not None
        #if isinstance(nodeID, string_types):
            #nodeID = 0
        self.element_layer[self.itotal, :] = [eid, layer]
        self.data[self.itime, self.itotal, :] = [o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.nelements
        ntimes = self.ntimes
        #nnodes = self.nnodes
        ntotal = self.ntotal
        nelements = len(unique(self.element_layer[:, 0]))

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i ntotal=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, ntotal))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i ntotal=%i\n'
                       % (self.__class__.__name__, nelements, ntotal))
            ntimes_word = 1
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        ctria3_msg, ctria6_msg, cquad4_msg, cquad8_msg = self._get_msgs()

        if self.element_type == 95:  # CQUAD4
            msg = cquad4_msg
            nnodes = 4
        elif self.element_type == 96:  # CQUAD8
            msg = cquad8_msg
            nnodes = 8
        elif self.element_type == 97:  # CTRIA3
            msg = ctria3_msg
            nnodes = 3
        elif self.element_type == 98:  # CTRIA6
            msg = ctria6_msg
            nnodes = 6
        else:
            raise NotImplementedError(self.element_name)

        return self.element_name, nnodes, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_layer[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element_layer[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #msg, nnodes, is_bilinear = self._get_msgs()
        if self.is_von_mises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        if self.is_strain():
            words = ['   ELEMENT  PLY   STRAINS IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR   STRAINS  PRINCIPAL  STRAINS (ZERO SHEAR)      %s\n' % von,
                     '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]
        else:
            words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' % von,
                     '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]

        if self.element_type == 95:  # CQUAD4
            if self.is_strain():
                msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
            else:
                msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
        #elif self.element_type == 96:  # CQUAD8
            #nnodes_per_element = 1
        elif self.element_type == 97:  # CTRIA3
            if self.is_strain():
                msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
            else:
                msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        elif self.element_type == 96:  # QUAD8
            # good
            if self.is_strain():
                msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words
            else:
                msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words

        elif self.element_type == 98:  # CTRIA6
            # good
            if self.is_strain():
                msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
            else:
                msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_layer[:, 0]
        layers = self.element_layer[:, 1]

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
            o11 = self.data[itime, :, 0]
            o22 = self.data[itime, :, 1]
            t12 = self.data[itime, :, 2]
            t1z = self.data[itime, :, 3]
            t2z = self.data[itime, :, 4]
            angle = self.data[itime, :, 5]
            major = self.data[itime, :, 6]
            minor = self.data[itime, :, 7]
            ovm = self.data[itime, :, 8]

            # loop over all the elements
            for eid, layer, o11i, o22i, t12i, t1zi, t2zi, anglei, majori, minori, ovmi in zip(
                eids, layers, o11, o22, t12, t1z, t2z, angle, major, minor, ovm):

                ([o11i, o22i, t12i, t1zi, t2zi, majori, minori, ovmi], is_all_zeros) = writeFloats12E([
                  o11i, o22i, t12i, t1zi, t2zi, majori, minori, ovmi])
                f.write('0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %s\n'
                        % (eid, layer, o11i, o22i, t12i, t1zi, t2zi, anglei, majori, minori, ovmi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealCompositePlateStressArray(RealCompositePlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCompositePlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def isStress(self):
        return True

    def isStrain(self):
        return False

    def get_headers(self):
        if self.is_von_mises():
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = ['o11', 'o22', 't12', 't1z', 't2z', 'angle', 'major', 'minor', ovm]
        return headers


class RealCompositePlateStrainArray(RealCompositePlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCompositePlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def isStress(self):
        return False

    def isStrain(self):
        return True

    def get_headers(self):
        if self.is_von_mises():
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = ['e11', 'e22', 'e12', 'e1z', 'e2z', 'angle', 'major', 'minor', ovm]
        return headers


class RealCompositePlateStress(StressObject):
    """
    ::

      # s_code = 0
                      S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
      ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
        ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.o11 = {}
        self.o22 = {}
        self.t12 = {}
        self.t1z = {}
        self.t2z = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        #self.fiberCurvature = {}
        self.ovmShear = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None, dt
            raise NotImplementedError('SORT2')
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.o11)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, fiberCurvature, o11, o22, t12, t1z, t2z, angle, '
                   'majorP, minorP, ovmShear\n')
        return msg

    def add_f06_data(self, data, transient, eType):
        """
        ::

                         S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
          ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
            ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
              151    1  -1.02406E+04  4.18348E+05  4.14359E+02   -8.62021E+00  1.86352E+04   89.94  4.18348E+05 -1.02410E+04  2.14295E+05
        """
        if transient is None:
            for line in data:
                if eType in ['CQUAD4', 'CTRIA3']:
                    (eid, ilayer, o11, o22, t12, t1z, t2z, angle, majorP, minorP, ovmShear) = line
                    if eid not in self.eType:
                        self.eType[eid] = eType
                        self.o11[eid] = [o11]
                        self.o22[eid] = [o22]
                        self.t12[eid] = [t12]
                        self.t1z[eid] = [t1z]
                        self.t2z[eid] = [t2z]
                        self.angle[eid] = [angle]
                        self.majorP[eid] = [majorP]
                        self.minorP[eid] = [minorP]
                        self.ovmShear[eid] = [ovmShear]
                    else:
                        self.o11[eid].append(o11)
                        self.o22[eid].append(o22)
                        self.t12[eid].append(t12)
                        self.t1z[eid].append(t1z)
                        self.t2z[eid].append(t2z)
                        self.angle[eid].append(angle)
                        self.majorP[eid].append(majorP)
                        self.minorP[eid].append(minorP)
                        self.ovmShear[eid].append(ovmShear)
                else:
                    msg = 'etype=%r is not supported...' % eType
                    raise NotImplementedError(msg)
            return
        raise NotImplementedError('transient results not supported')

    def delete_transient(self, dt):
        del self.o11[dt]
        del self.o22[dt]
        del self.t12[dt]
        del self.t1z[dt]
        del self.t2z[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.ovmShear[dt]

    def get_transients(self):
        k = self.o11.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #self.fiberDistance[dt] = {}
        self.o11[dt] = {}
        self.o22[dt] = {}
        self.t12[dt] = {}
        self.t1z[dt] = {}
        self.t2z[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.ovmShear[dt] = {}

    def add_new_eid(self, eType, dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                    majorP, minorP, ovm):
        """all points are located at the centroid"""
        if eid in self.o11:
            return self.add(dt, eid, layer, o11, o22, t12, t1z, t2z, angle, majorP, minorP, ovm)
        assert eid not in self.o11, 'eid=%s o11=%s eType=%s code=%s' % (
            eid, self.o11[eid], self.eType[eid], self.data_code)

        self.eType[eid] = eType
        self.o11[eid] = [o11]
        self.o22[eid] = [o22]
        self.t12[eid] = [t12]
        self.t1z[eid] = [t1z]
        self.t2z[eid] = [t2z]
        self.angle[eid] = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.ovmShear[eid] = [ovm]
        #if nodeID==0: raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                          majorP, minorP, ovm):
        """all points are located at the centroid"""
        if dt not in self.o11:
            self.add_new_transient(dt)
        if eid in self.o11[dt]:
            return self.add(dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                            majorP, minorP, ovm)

        assert eid not in self.o11[dt]
        self.eType[eid] = eType
        self.o11[dt][eid] = [o11]
        self.o22[dt][eid] = [o22]
        self.t12[dt][eid] = [t12]
        self.t1z[dt][eid] = [t1z]
        self.t2z[dt][eid] = [t2z]
        self.angle[dt][eid] = [angle]
        self.majorP[dt][eid] = [majorP]
        self.minorP[dt][eid] = [minorP]
        self.ovmShear[dt][eid] = [ovm]
        #if nodeID==0: raise Exception(msg)

    def add(self, dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
            majorP, minorP, ovm):
        self.o11[eid].append(o11)
        self.o22[eid].append(o22)
        self.t12[eid].append(t12)
        self.t1z[eid].append(t1z)
        self.t2z[eid].append(t2z)
        self.angle[eid].append(angle)
        self.majorP[eid].append(majorP)
        self.minorP[eid].append(minorP)
        self.ovmShear[eid].append(ovm)
        #if nodeID==0: raise Exception(msg)

    def add_sort1(self, dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                 majorP, minorP, ovm):
        self.o11[dt][eid].append(o11)
        self.o22[dt][eid].append(o22)
        self.t12[dt][eid].append(t12)
        self.t1z[dt][eid].append(t1z)
        self.t2z[dt][eid].append(t2z)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)
        self.ovmShear[dt][eid].append(ovm)

    def getHeaders(self):
        headers = ['o11', 'o22', 't12', 't1z', 't2z']
        if self.is_von_mises():
            headers.append('oVonMises')
        else:
            headers.append('maxShear')
        return headers

    def _is_cquad4(self, etype):
        if self.element_type == 95:
        #if etype in ['CQUAD4', 'QUAD4LC', 'QUAD4LC-composite']:
            return True
        elif self.element_type == 97:
        #elif etype in ['CTRIA3', 'TRIA3LC', 'TRIA3LC-composite']:
            return False
        elif self.element_type == 96:
        #elif etype in ['QUAD8LC-composite']:
            return False
        elif self.element_type == 98:
        #elif etype in ['TRIA6LC-composite']:
            return False
        msg = 'element_type=%s element_number=%s' % (etype, self.element_type)
        raise NotImplementedError(msg)

    def _is_ctria3(self, etype):
        if self.element_type == 95:
        #if etype in ['CQUAD4', 'QUAD4LC', 'QUAD4LC-composite']:
            return False
        elif self.element_type == 97:
        #elif etype in ['CTRIA3', 'TRIA3LC', 'TRIA3LC-composite']:
            return True
        elif self.element_type == 96:
        #elif etype in ['QUAD8LC-composite']:
            return False
        elif self.element_type == 98:
        #elif etype in ['TRIA6LC-composite']:
            return False
        msg = 'element_type=%s element_number=%s' % (etype, self.element_type)
        raise NotImplementedError(msg)

    def _is_cquad8(self, etype):
        if self.element_type == 95:
        #if etype in ['CQUAD4', 'QUAD4LC', 'QUAD4LC-composite']:
            return False
        elif self.element_type == 97:
        #elif etype in ['CTRIA3', 'TRIA3LC', 'TRIA3LC-composite']:
            return False
        elif self.element_type == 96:
        #elif etype in ['QUAD8LC-composite']:
            return True
        elif self.element_type == 98:
        #elif etype in ['TRIA6LC-composite']:
            return False
        msg = 'element_type=%s element_number=%s' % (etype, self.element_type)
        raise NotImplementedError(msg)

    def _is_ctria6(self, etype):
        if self.element_type == 95:  # 33, 144,
        #if etype in ['CQUAD4', 'QUAD4LC', 'QUAD4LC-composite']:
            return False
        elif self.element_type == 97: # 74,
        #elif etype in ['CTRIA3', 'TRIA3LC', 'TRIA3LC-composite']:
            return False
        elif self.element_type == 96:
        #elif etype in ['QUAD8LC-composite']:
            return False
        elif self.element_type == 98:
        #elif etype in ['TRIA6LC-composite']:
            return True
        msg = 'element_type=%s element_number=%s' % (etype, self.element_type)
        raise NotImplementedError(msg)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        if self.is_von_mises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' % von,
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]

        etypes = list(self.eType.values())
        if self._is_cquad4(etypes[0]):
            cquad4_msg = header + ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
        elif self._is_ctria3(etypes[0]):
            cquad4_msg = header + ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        elif self._is_ctria6(etypes[0]):
            cquad4_msg = header + ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
        elif self._is_cquad8(etypes[0]):
            cquad4_msg = header + ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words
        else:
            raise NotImplementedError(etypes)

        for eid, o11s in sorted(iteritems(self.o11)):
            out = ''
            etype = self.eType[eid]
            for ilayer in range(len(o11s)):
                o11 = self.o11[eid][ilayer]
                o22 = self.o22[eid][ilayer]
                t12 = self.t12[eid][ilayer]
                t1z = self.t1z[eid][ilayer]
                t2z = self.t2z[eid][ilayer]

                angle = self.angle[eid][ilayer]
                major = self.majorP[eid][ilayer]
                minor = self.minorP[eid][ilayer]
                ovm = self.ovmShear[eid][ilayer]
                (vals2, is_all_zeros) = writeFloats12E([o11, o22, t12, t1z, t2z, major, minor, ovm])
                [o11, o22, t12, t1z, t2z, major, minor, ovm] = vals2
                out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' % (eid, ilayer + 1, o11, o22, t12, t1z, t2z, angle, major, minor, ovm)

            cquad4_msg.append(out)

        cquad4_msg.append(page_stamp % page_num)
        page_num += 1

        f.write(''.join(cquad4_msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             is_mag_phase=False, is_sort1=True):
        if self.is_von_mises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' % von,
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]

        etypes = list(self.eType.values())
        if self._is_cquad4(list(etypes)[0]):
            cquad4_words = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
        elif self._is_ctria3(list(etypes)[0]):
            cquad4_words = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        elif self._is_ctria6(list(etypes)[0]):
            cquad4_words = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
        elif self._is_cquad8(list(etypes)[0]):
            cquad4_words = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words
        else:
            raise NotImplementedError(etypes)

        for dt, O11s in sorted(iteritems(self.o11)):
            quad_msg = []
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            quad_msg = header + cquad4_words

            for eid, o11s in sorted(iteritems(O11s)):
                out = ''
                etype = self.eType[eid]
                for ilayer in range(len(o11s)):
                    o11 = self.o11[dt][eid][ilayer]
                    o22 = self.o22[dt][eid][ilayer]
                    t12 = self.t12[dt][eid][ilayer]
                    t1z = self.t1z[dt][eid][ilayer]
                    t2z = self.t2z[dt][eid][ilayer]

                    angle = self.angle[dt][eid][ilayer]
                    major = self.majorP[dt][eid][ilayer]
                    minor = self.minorP[dt][eid][ilayer]
                    ovm = self.ovmShear[dt][eid][ilayer]
                    (vals2, is_all_zeros) = writeFloats12E([o11, o22,
                                                          t12, t1z, t2z,
                                                          major, minor, ovm])
                    [o11, o22, t12, t1z, t2z, major, minor, ovm] = vals2
                    out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' % (eid, ilayer + 1, o11, o22, t12, t1z, t2z, angle, major, minor, ovm)
                quad_msg.append(out)

            quad_msg.append(page_stamp % page_num)
            f.write(''.join(quad_msg))
            page_num += 1
        return page_num - 1


class RealCompositePlateStrain(StrainObject):
    """
    ::

      ???
      ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
        ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)

        self.eType = {}
        self.code = [self.format_code, self.sort_code, self.s_code]
        self.e11 = {}
        self.e22 = {}
        self.e12 = {}
        self.e1z = {}
        self.e2z = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}

        if self.code == [1, 0, 14]:
            self.evmShear = {}
            assert self.is_von_mises() == False, 'format/sort/stressCode=%s bits=%s' % (self.code, self.stress_bits)
        else:
            raise RuntimeError('Invalid Code: compositePlateStrain - get the format/sort/stressCode=%s bits=%s' % (self.code, self.stress_bits))

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None, dt
            raise NotImplementedError('SORT2')
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.e11)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, e11, e22, e12, e1z, e2z, angle, majorP, minorP\n')
        return msg

    def add_f06_data(self, data, transient, eType):
        """
        ::

                         S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
          ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
            ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
              151    1  -1.02406E+04  4.18348E+05  4.14359E+02   -8.62021E+00  1.86352E+04   89.94  4.18348E+05 -1.02410E+04  2.14295E+05
        """
        if transient is None:
            for line in data:
                if eType in ['CQUAD4', 'CTRIA3']:
                    (eid, ilayer, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evmShear) = line
                    if eid not in self.eType:
                        self.eType[eid] = eType
                        self.e11[eid] = [e11]
                        self.e22[eid] = [e22]
                        self.e12[eid] = [e12]
                        self.e1z[eid] = [e1z]
                        self.e2z[eid] = [e2z]
                        self.angle[eid] = [angle]
                        self.majorP[eid] = [majorP]
                        self.minorP[eid] = [minorP]
                        self.evmShear[eid] = [evmShear]
                    else:
                        self.e11[eid].append(e11)
                        self.e22[eid].append(e22)
                        self.e12[eid].append(e12)
                        self.e1z[eid].append(e1z)
                        self.e2z[eid].append(e2z)
                        self.angle[eid].append(angle)
                        self.majorP[eid].append(majorP)
                        self.minorP[eid].append(minorP)
                        self.evmShear[eid].append(evmShear)
                else:
                    msg = 'etype=%r is not supported...' % eType
                    raise NotImplementedError(msg)
            return
        raise NotImplementedError('transient results not supported')

    def delete_transient(self, dt):
        del self.e11[dt]
        del self.e22[dt]
        del self.e12[dt]
        del self.e1z[dt]
        del self.e2z[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.evmShear[dt]

    def get_transients(self):
        k = self.e11.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #self.fiberDistance[dt] = {}
        self.e11[dt] = {}
        self.e22[dt] = {}
        self.e12[dt] = {}
        self.e1z[dt] = {}
        self.e2z[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.evmShear[dt] = {}

    def add_new_eid(self, eType, dt, eid, layer, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm):
        """all points are located at the centroid"""
        if eid in self.e11:
            return self.add(dt, eid, layer, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm)
        assert eid not in self.e11, eid
        assert isinstance(eid, int), eid
        self.eType[eid] = eType
        self.e11[eid] = [e11]
        self.e22[eid] = [e22]
        self.e12[eid] = [e12]
        self.e1z[eid] = [e1z]
        self.e2z[eid] = [e2z]
        self.angle[eid] = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.evmShear[eid] = [evm]
        #msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid, layer, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm):
        """all points are located at the centroid"""
        if dt not in self.e11:
            self.add_new_transient(dt)
        assert eid not in self.e11[dt], eid
        assert isinstance(eid, int), eid

        self.eType[eid] = eType
        self.e11[dt][eid] = [e11]
        self.e22[dt][eid] = [e22]
        self.e12[dt][eid] = [e12]
        self.e1z[dt][eid] = [e1z]
        self.e2z[dt][eid] = [e2z]
        self.angle[dt][eid] = [angle]
        self.majorP[dt][eid] = [majorP]
        self.minorP[dt][eid] = [minorP]
        self.evmShear[dt][eid] = [evm]

    def add(self, dt, eid, layer, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm):
        self.e11[eid].append(e11)
        self.e22[eid].append(e22)
        self.e12[eid].append(e12)
        self.e1z[eid].append(e1z)
        self.e2z[eid].append(e2z)
        self.angle[eid].append(angle)
        self.majorP[eid].append(majorP)
        self.minorP[eid].append(minorP)
        self.evmShear[eid].append(evm)
        #if nodeID==0: raise Exception(msg)

    def add_sort1(self, dt, eid, layer, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm):
        self.e11[dt][eid].append(e11)
        self.e22[dt][eid].append(e22)
        self.e12[dt][eid].append(e12)
        self.e1z[dt][eid].append(e1z)
        self.e2z[dt][eid].append(e2z)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)
        self.evmShear[dt][eid].append(evm)
        #if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        headers = ['e11', 'e22', 'e12', 'e1z', 'e2z']
        if self.is_von_mises():
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        if self.is_von_mises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY   STRAINS IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR   STRAINS  PRINCIPAL  STRAINS (ZERO SHEAR)      %s\n' % von,
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]

        etypes = list(self.eType.values())
        if self._is_cquad4(list(etypes)[0]):
            cquad4_msg = header + ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
        elif self._is_ctria3(list(etypes)[0]):
            cquad4_msg = header + ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        elif self._is_ctria6(list(etypes)[0]):
            cquad4_msg = header + ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
        elif self._is_cquad8(list(etypes)[0]):
            cquad4_msg = header + ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words
        else:
            raise NotImplementedError(etypes)

        for eid, e11s in sorted(iteritems(self.e11)):
            out = ''
            #etype = self.eType[eid]
            for ilayer in range(len(e11s)):
                e11 = self.e11[eid][ilayer]
                e22 = self.e22[eid][ilayer]
                e12 = self.e12[eid][ilayer]
                e1z = self.e1z[eid][ilayer]
                e2z = self.e2z[eid][ilayer]

                angle = self.angle[eid][ilayer]
                major = self.majorP[eid][ilayer]
                minor = self.minorP[eid][ilayer]
                evm = self.evmShear[eid][ilayer]

                (vals2, is_all_zeros) = writeFloats12E([e11,
                                                        e22, e12, e1z, e2z, major, minor, evm])
                [e11, e22, e12, e1z, e2z, major, minor, evm] = vals2
                out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' % (eid, ilayer + 1, e11, e22, e12, e1z, e2z, angle, major, minor, evm)
            cquad4_msg.append(out)
        cquad4_msg.append(page_stamp % page_num)
        page_num += 1

        f.write(''.join(cquad4_msg))
        return page_num

    def _is_cquad4(self, etype):
        if self.element_type == 95:
        #if etype in ['CQUAD4', 'QUAD4LC', 'QUAD4LC-composite']:
            return True
        elif self.element_type == 97:
        #elif etype in ['CTRIA3', 'TRIA3LC', 'TRIA3LC-composite']:
            return False
        elif self.element_type == 96:
        #elif etype in ['QUAD8LC-composite']:
            return False
        elif self.element_type == 98:
        #elif etype in ['TRIA6LC-composite']:
            return False
        msg = 'element_type=%s element_number=%s' % (etype, self.element_type)
        raise NotImplementedError(msg)

    def _is_ctria3(self, etype):
        if self.element_type == 95:
        #if etype in ['CQUAD4', 'QUAD4LC', 'QUAD4LC-composite']:
            return False
        elif self.element_type == 97:
        #elif etype in ['CTRIA3', 'TRIA3LC', 'TRIA3LC-composite']:
            return True
        elif self.element_type == 96:
        #elif etype in ['QUAD8LC-composite']:
            return False
        elif self.element_type == 98:
        #elif etype in ['TRIA6LC-composite']:
            return False
        msg = 'element_type=%s element_number=%s' % (etype, self.element_type)
        raise NotImplementedError(msg)

    def _is_cquad8(self, etype):
        if self.element_type == 95:
        #if etype in ['CQUAD4', 'QUAD4LC', 'QUAD4LC-composite']:
            return False
        elif self.element_type == 97:
        #elif etype in ['CTRIA3', 'TRIA3LC', 'TRIA3LC-composite']:
            return False
        elif self.element_type == 96:
        #elif etype in ['QUAD8LC-composite']:
            return True
        elif self.element_type == 98:
        #elif etype in ['TRIA6LC-composite']:
            return False
        msg = 'element_type=%s element_number=%s' % (etype, self.element_type)
        raise NotImplementedError(msg)

    def _is_ctria6(self, etype):
        if self.element_type == 95:
        #if etype in ['CQUAD4', 'QUAD4LC', 'QUAD4LC-composite']:
            return False
        elif self.element_type == 97:
            return False
        #elif etype in ['CTRIA3', 'TRIA3LC', 'TRIA3LC-composite']:
        elif self.element_type == 96:
        #elif etype in ['QUAD8LC-composite']:
            return False
        elif self.element_type == 98:
        #elif etype in ['TRIA6LC-composite']:
            return True
        msg = 'element_type=%s element_number=%s' % (etype, self.element_type)
        raise NotImplementedError(msg)

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             is_mag_phase=False, is_sort1=True):
        if self.is_von_mises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        common_words = ['   ELEMENT  PLY   STRAINS IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR   STRAINS  PRINCIPAL  STRAINS (ZERO SHEAR)      %s\n' % von,
                        '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]
        words = []

        etypes = set(self.eType.values())
        if self._is_cquad4(list(etypes)[0]):
            words = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + common_words
        elif self._is_ctria3(list(etypes)[0]):
            words = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + common_words
        elif self._is_ctria6(list(etypes)[0]):
            words = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + common_words
        elif self._is_cquad8(list(etypes)[0]):
            words = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + common_words
        else:
            raise NotImplementedError(etypes)

        for dt, e11s in sorted(iteritems(self.e11)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            quad_msg = header + words
            for eid, e11s in sorted(iteritems(e11s)):
                out = ''
                #etype = self.eType[eid]
                for iLayer in range(len(e11s)):
                    e11 = self.e11[dt][eid][iLayer]
                    e22 = self.e22[dt][eid][iLayer]
                    e12 = self.e12[dt][eid][iLayer]
                    e1z = self.e1z[dt][eid][iLayer]
                    e2z = self.e2z[dt][eid][iLayer]

                    angle = self.angle[dt][eid][iLayer]
                    major = self.majorP[dt][eid][iLayer]
                    minor = self.minorP[dt][eid][iLayer]
                    evm = self.evmShear[dt][eid][iLayer]
                    (vals2, is_all_zeros) = writeFloats12E([e11, e22,
                                                            e12, e1z, e2z, major, minor, evm])
                    [e11, e22, e12, e1z, e2z, major, minor, evm] = vals2
                    out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' % (eid, iLayer + 1, e11, e22, e12, e1z, e2z, angle, major, minor, evm)
                quad_msg.append(out)

            quad_msg.append(page_stamp % page_num)
            page_num += 1
            f.write(''.join(quad_msg))
        return page_num - 1
