from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from itertools import count
from numpy import zeros, searchsorted, ravel

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, _eigenvalue_header, get_key0


class RealTriaxArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific

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
        return headers

    def build(self):
        if self.is_built:
            return
        #print("self.ielement =", self.ielement)
        print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        #aasdf
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        if self.element_type == 53:
            nnodes_per_element = 1
        else:
            raise NotImplementedError(self.element_type)

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
        print('self.element_node.shape', self.element_node.shape)

        # [radial, azimuthal, axial, shear, omax, oms, ovm]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

    def add_sort1(self, dt, eid, nid, radial, azimuthal, axial, shear, omax, oms, ovm):
        assert isinstance(eid, int)
        self._times[self.itime] = dt
        self.element_node[self.itotal, :] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [radial, azimuthal, axial, shear, omax, oms, ovm]
        self.itotal += 1
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.ntotal
        ntimes = self.ntimes
        ntotal = self.ntotal
        nelements = self.ntotal

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element == eid) for eid in eids])
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = self._get_msgs()
        (ntimes, ntotal) = self.data.shape[:2]
        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            #[radial, azimuthal, axial, shear, omax, oms, ovm]
            radial = self.data[itime, :, 0]
            azimuthal = self.data[itime, :, 1]
            axial = self.data[itime, :, 2]
            shear = self.data[itime, :, 3]
            omax = self.data[itime, :, 4]
            oms = self.data[itime, :, 5]
            ovm = self.data[itime, :, 6]

            # loop over all the elements
            for (i, eid, nid, radiali, azimuthali, axiali, sheari, omaxi, omsi, ovmi) in zip(
                count(), eids, nids, radial, azimuthal, axial, shear, omax, oms, ovm):

                vals = [radiali, azimuthali, axiali, sheari, omaxi, omsi, ovmi]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [radiali, azimuthali, axiali, sheari, omaxi, omsi, ovmi] = vals2
                f.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                    % (eid, nid, radiali, azimuthali, axiali, sheari, omaxi, omsi, ovmi))
            f.write(page_stamp % page_num)
            page_num += 1
        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num


class RealTriaxStressArray(RealTriaxArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTriaxArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['radial', 'azimuthal', 'axial', 'shear', 'omax', 'oms', 'ovm']
        return headers

    def _get_msgs(self):
        if self.element_type == 53:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = ['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
               '   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
               '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n',
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
        ]
        return msg

class RealTriaxStrainArray(RealTriaxArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTriaxArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['radial', 'azimuthal', 'axial', 'shear', 'omax', 'oms', 'ovm']
        return headers

    def _get_msgs(self):
        if self.element_type == 53:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = ['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
               '   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
               '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n',
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
        ]
        return msg

#class RealTriaxStress(StressObject):
    #"""
    #::

      ## format_code=1 sort_code=0 stressCode=0
                                        #S T R E S S E S   I N   T R I A X 6   E L E M E N T S
      #ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES
         #ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR
         #5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                  #4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
    #"""
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StressObject.__init__(self, data_code, isubcase)
        #self.eType = 'CTRIAX6'

        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.radial = {}
        #self.azimuthal = {}
        #self.axial = {}
        #self.shear = {}
        #self.omax = {}
        #self.oms = {}
        #self.ovm = {}

        #self.dt = dt
        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        #else:
            #assert dt is not None
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2

    #def is_real(self):
        #return True

    #def is_complex(self):
        #return False

    #def get_stats(self):
        #msg = self.get_data_code()
        #if self.nonlinear_factor is not None:  # transient
            #ntimes = len(self.radial)
            #r0 = get_key0(self.radial)
            #nelements = len(self.radial[r0])
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #nelements = len(self.radial)
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
        #msg.append('  eType, radial, azimuthal, axial, shear, '
                   #'omax, oms, ovm\n')
        #return msg

    #def add_f06_data(self, data, transient):
        #raise Exception('Not Implemented')
        #if transient is None:
            #for line in data:
                #(eid, axial, MSa, torsion, MSt) = line
                #if MSa is None:
                    #MSa = 0.
                #if MSt is None:
                    #MSt = 0.
                #self.axial[eid] = axial
                #self.MS_axial[eid] = MSa
                #self.torsion[eid] = torsion
                #self.MS_torsion[eid] = MSt
            #return

        #(dtName, dt) = transient
        #self.data_code['name'] = dtName
        #if dt not in self.s1:
            #self.update_dt(self.data_code, dt)
            #self.isTransient = True

        #for line in data:
            #(eid, axial, MSa, torsion, MSt) = line
            #if MSa is None:
                #MSa = 0.
            #if MSt is None:
                #MSt = 0.
            #self.axial[dt][eid] = axial
            #self.MS_axial[dt][eid] = MSa
            #self.torsion[dt][eid] = torsion
            #self.MS_torsion[dt][eid] = MSt

    #def delete_transient(self, dt):
        #del self.radial[dt]
        #del self.azimuthal[dt]
        #del self.axial[dt]
        #del self.shear[dt]
        #del self.omax[dt]
        #del self.oms[dt]
        #del self.ovm[dt]

    #def get_transients(self):
        #k = self.axial.keys()
        #k.sort()
        #return k

    #def add_new_transient(self, dt):
        #"""
        #initializes the transient variables
        #"""
        #self.radial[dt] = {}
        #self.azimuthal[dt] = {}
        #self.axial[dt] = {}
        #self.shear[dt] = {}
        #self.omax[dt] = {}
        #self.oms[dt] = {}
        #self.ovm[dt] = {}

    #def add_new_eid(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        ##print "**?eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        #self.radial[eid] = {nid: rs}
        #self.azimuthal[eid] = {nid: azs}
        #self.axial[eid] = {nid: As}
        #self.shear[eid] = {nid: ss}
        #self.omax[eid] = {nid: maxp}
        #self.oms[eid] = {nid: tmax}
        #self.ovm[eid] = {nid: octs}

    #def add(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        ##print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        #self.radial[eid][nid] = rs
        #self.azimuthal[eid][nid] = azs
        #self.axial[eid][nid] = As
        #self.shear[eid][nid] = ss
        #self.omax[eid][nid] = maxp
        #self.oms[eid][nid] = tmax
        #self.ovm[eid][nid] = octs

    #def add_new_eid_sort1(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        ##assert isinstance(eid,int)
        ##assert eid >= 0, eid
        ##print "*  eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        #if dt not in self.radial:
            #self.add_new_transient(dt)
        #self.radial[dt][eid] = {nid: rs}
        #self.azimuthal[dt][eid] = {nid: azs}
        #self.axial[dt][eid] = {nid: As}
        #self.shear[dt][eid] = {nid: ss}
        #self.omax[dt][eid] = {nid: maxp}
        #self.oms[dt][eid] = {nid: tmax}
        #self.ovm[dt][eid] = {nid: octs}

    #def add_sort1(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        ##print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        #self.radial[dt][eid][nid] = rs
        #self.azimuthal[dt][eid][nid] = azs
        #self.axial[dt][eid][nid] = As
        #self.shear[dt][eid][nid] = ss
        #self.omax[dt][eid][nid] = maxp
        #self.oms[dt][eid][nid] = tmax
        #self.ovm[dt][eid][nid] = octs

    #def write_f06(self, header, page_stamp, page_num=1, f=None,
                  #is_mag_phase=False, is_sort1=True):
        #if self.nonlinear_factor is not None:
            #return self._write_f06_transient(header, page_stamp, page_num, f,
                                             #is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        #msg = header + ['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
                        #'   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                        #'      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
                       ##'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                       ##'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        ##out = []
        #for eid, radial in sorted(iteritems(self.radial)):
            #for nid in sorted(radial):
                #rad = self.radial[eid][nid]
                #azimuth = self.azimuthal[eid][nid]
                #axial = self.axial[eid][nid]
                #shear = self.shear[eid][nid]
                #omax = self.omax[eid][nid]
                #oms = self.oms[eid][nid]
                #ovm = self.ovm[eid][nid]
                #if nid == 0:
                    #Eid = eid
                #else:
                    #Eid = ''
                #([rad, azimuth, axial, shear, omax, oms, ovm], is_all_zeros) = writeFloats13E([rad, azimuth, axial, shear, omax, oms, ovm])
                #msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, radial, azimuth, axial, shear, omax, oms, ovm.rstrip()))
            #msg.append('\n')

        #msg.append(page_stamp % page_num)
        #f.write(''.join(msg))
        #return page_num

    #def _write_f06_transient(self, header, page_stamp,
                             #page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #words = ['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
                 #'   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                 #'      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
                ##'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                ##'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        #msg = []
        #for dt, Radial in sorted(iteritems(self.radial)):
            #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            #msg += header + words
            #for eid, radial in sorted(iteritems(Radial)):
                #for nid in sorted(radial):
                    #rad = self.radial[dt][eid][nid]
                    #azimuth = self.azimuthal[dt][eid][nid]
                    #axial = self.axial[dt][eid][nid]
                    #shear = self.shear[dt][eid][nid]
                    #omax = self.omax[dt][eid][nid]
                    #oms = self.oms[dt][eid][nid]
                    #ovm = self.ovm[dt][eid][nid]
                    #if nid == 0:
                        #Eid = eid
                    #else:
                        #Eid = ''
                    #([rad, azimuth, axial, shear, omax, oms, ovm], is_all_zeros) = writeFloats13E([rad, azimuth, axial, shear, omax, oms, ovm])
                    #msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, rad, azimuth, axial, shear, omax, oms, ovm.rstrip()))
                #msg.append('\n')

            #msg.append(page_stamp % page_num)
            #f.write(''.join(msg))
            #msg = ['']
            #page_num += 1
        #return page_num - 1


#class RealTriaxStrain(StrainObject):
    #"""
    #::

      ## format_code=1 sort_code=0 stressCode=0
                                        #S T R A I N S   I N   T R I A X 6   E L E M E N T S
      #ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES
         #ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR
         #5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                  #4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
    #"""
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StrainObject.__init__(self, data_code, isubcase)
        #self.eType = 'CTRIAX6'

        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.radial = {}
        #self.azimuthal = {}
        #self.axial = {}
        #self.shear = {}
        #self.emax = {}
        #self.ems = {}
        #self.evm = {}

        #self.dt = dt
        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        #else:
            #assert dt is not None
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2

    #def get_stats(self):
        #msg = self.get_data_code()
        #if self.nonlinear_factor is not None:  # transient
            #ntimes = len(self.radial)
            #r0 = get_key0(self.radial)
            #nelements = len(self.radial[r0])
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #nelements = len(self.radial)
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
        #msg.append('  eType, radial, azimuthal, axial, shear, '
                   #'emax, ems, evm\n')
        #return msg

    #def add_f06_data(self, data, transient):
        #raise Exception('Not Implemented')

    #def delete_transient(self, dt):
        #del self.radial[dt]
        #del self.azimuthal[dt]
        #del self.axial[dt]
        #del self.shear[dt]
        #del self.emax[dt]
        #del self.ems[dt]
        #del self.evm[dt]

    #def get_transients(self):
        #k = self.axial.keys()
        #k.sort()
        #return k

    #def add_new_transient(self, dt):
        #"""
        #initializes the transient variables
        #"""
        #self.radial[dt] = {}
        #self.azimuthal[dt] = {}
        #self.axial[dt] = {}
        #self.shear[dt] = {}
        #self.emax[dt] = {}
        #self.ems[dt] = {}
        #self.evm[dt] = {}

    #def add_new_eid(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        ##print "**?eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        #self.radial[eid] = {nid: rs}
        #self.azimuthal[eid] = {nid: azs}
        #self.axial[eid] = {nid: As}
        #self.shear[eid] = {nid: ss}
        #self.emax[eid] = {nid: maxp}
        #self.ems[eid] = {nid: emax}
        #self.evm[eid] = {nid: ects}

    #def add(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        ##print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        #self.radial[eid][nid] = rs
        #self.azimuthal[eid][nid] = azs
        #self.axial[eid][nid] = As
        #self.shear[eid][nid] = ss
        #self.emax[eid][nid] = maxp
        #self.ems[eid][nid] = emax
        #self.evm[eid][nid] = ects

    #def add_new_eid_sort1(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        ##assert isinstance(eid,int)
        ##assert eid >= 0, eid
        ##print "*  eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        #self.radial[dt][eid] = {nid: rs}
        #self.azimuthal[dt][eid] = {nid: azs}
        #self.axial[dt][eid] = {nid: As}
        #self.shear[dt][eid] = {nid: ss}
        #self.emax[dt][eid] = {nid: maxp}
        #self.ems[dt][eid] = {nid: emax}
        #self.evm[dt][eid] = {nid: ects}

    #def add_sort1(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        ##print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        #self.radial[dt][eid][nid] = rs
        #self.azimuthal[dt][eid][nid] = azs
        #self.axial[dt][eid][nid] = As
        #self.shear[dt][eid][nid] = ss
        #self.emax[dt][eid][nid] = maxp
        #self.ems[dt][eid][nid] = emax
        #self.evm[dt][eid][nid] = ects

    #def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #if self.nonlinear_factor is not None:
            #return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        #msg = header + ['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
                        #'   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                        #'      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              ##         '      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              ##         '               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        ##out = []
        #for eid, radial in sorted(iteritems(self.radial)):
            #for nid in sorted(radial):
                #rad = self.radial[eid][nid]
                #azimuth = self.azimuthal[eid][nid]
                #axial = self.axial[eid][nid]
                #shear = self.shear[eid][nid]
                #emax = self.emax[eid][nid]
                #ems = self.ems[eid][nid]
                #evm = self.evm[eid][nid]
                #if nid == 0:
                    #Eid = eid
                #else:
                    #Eid = ''
                #([rad, azimuth, axial, shear, emax, ems, evm], is_all_zeros) = writeFloats13E([rad, azimuth, axial, shear, emax, ems, evm])
                #msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, radial, azimuth, axial, shear, emax, ems, evm.rstrip()))
            #msg.append('\n')

        #msg.append(page_stamp % page_num)
        #f.write(''.join(msg))
        #return page_num

    #def _write_f06_transient(self, header, page_stamp,
                          #page_num=1, f=None, is_mag_phase=False):
        #words = ['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
                 #'   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                 #'      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              ##'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              ##'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        #msg = []
        #for dt, Radial in sorted(iteritems(self.radial)):
            #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            #msg += header + words
            #for eid, radial in sorted(iteritems(Radial)):
                #for nid in sorted(radial):
                    #rad = self.radial[dt][eid][nid]
                    #azimuth = self.azimuthal[dt][eid][nid]
                    #axial = self.axial[dt][eid][nid]
                    #shear = self.shear[dt][eid][nid]
                    #emax = self.emax[dt][eid][nid]
                    #ems = self.ems[dt][eid][nid]
                    #evm = self.evm[dt][eid][nid]
                    #if nid == 0:
                        #Eid = eid
                    #else:
                        #Eid = ''
                    #([rad, azimuth, axial, shear, emax, ems, evm], is_all_zeros) = writeFloats13E([rad, azimuth, axial, shear, emax, ems, evm])
                    #msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n'
                               #% (Eid, nid, rad, azimuth, axial, shear, emax,
                                  #ems, evm.rstrip()))
                #msg.append('\n')

            #msg.append(page_stamp % page_num)
            #f.write(''.join(msg))
            #msg = ['']
            #page_num += 1
        #return page_num - 1
