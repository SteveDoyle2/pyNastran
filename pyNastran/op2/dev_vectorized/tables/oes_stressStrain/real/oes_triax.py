from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import zeros, array
import pandas as pd

from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E

class RealTriaxObject(object):
    def __init__(self):
        self.shape = {}
        self._inode_start = None
        self._inode_end = None
        self._ielement_start = None
        self._ielement_end = None
        self._ncount = 0
        self.data = None
        self.element_data = None

    def _increase_size(self, dt, nnodes, nelements):
        if dt in self.shape:  # default dictionary
            self.shape[dt][0] += nnodes
            self.shape[dt][1] += nelements
        else:
            self.shape[dt] = [nnodes, nelements]

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nnodes = self.shape[shape0][0]
        nelements = self.shape[shape0][1]
        #print("ndt=%s nnodes=%s dts=%s" % (str(ndt), str(nnodes), str(dts)))
        return ndt, nnodes, nelements, dts

    def _preallocate(self, dt, nnodes, nelements):
        """
        nodes will create self.data
        elements will create self.element_data
        
        data is primary, so nodes are input first
        """
        #print('---preallocate nnodes=%s nelements=%s' % (nnodes, nelements))
        if self.shape is None:
            self._inode_start += nnodes
            self._inode_end += nnodes
            self._ielement_start += nelements
            self._ielement_end += nelements
        else:
            ndt, nnodes_size, nelements_size, dts = self._get_shape()
            #print("ndt=%s nelements_size=%s nnodes_size=%s dts=%s" % (ndt, nelements_size, nnodes_size, str(dts)))

            if self._inode_start is not None:
                return (self._inode_start, self._inode_start + nnodes,
                        self._ielement_start, self._ielement_start + nelements)

            n = ndt * nnodes_size
            if self._ncount != 0:
                asfd
            self._ncount += 1
            self._inode_start = 0
            self._inode_end = nnodes

            self._ielement_start = 0
            self._ielement_end = nelements

            data = {}
            #element_data = {}
            columns = []
            if dts[0] is not None:
                name = self.data_code['name']
                #print('***name=%r***' % name)
                if isinstance(dt, int):
                    data[name] = pd.Series(zeros((n), dtype='int32'))
                else:
                    data[name] = pd.Series(zeros((n), dtype='float32'))
                columns.append(name)

            #element_data['element_id'] = pd.Series(zeros((nelements_size), dtype='int32'))
            #element_data['element_type'] = pd.Series(zeros(nelements_size, dtype='str'))
            #element_data['nlayers'] = pd.Series(zeros(nelements_size, dtype='int32'))

            headers = self._get_headers()
            #(radial, azimuthal, axial, shear, omax, oms, evm) = headers

            data['element_id'] = pd.Series(zeros((n), dtype='int32'))
            data['node_id']   = pd.Series(zeros((n), dtype='int32'))
            data['azs']  = pd.Series(zeros((n), dtype='float32'))
            data['As']   = pd.Series(zeros((n), dtype='float32'))
            data['ss']   = pd.Series(zeros((n), dtype='float32'))
            data['maxp'] = pd.Series(zeros((n), dtype='float32'))
            data['tmax'] = pd.Series(zeros((n), dtype='float32'))
            data['octs'] = pd.Series(zeros((n), dtype='float32'))

            #columns.append('element_type')

            #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
            #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
            #print('n =', n)


            #headers = self._get_headers()
            #(fd, oxx, oyy, txy, omax, omin, ovm) = headers
            columns += ['element_id', 'node_id', 'azs', 'As', 'ss', 'maxp', 'tmax', 'octs']

            self.data = pd.DataFrame(data, columns=columns)
            #self.element_data = pd.DataFrame(element_data, columns=['element_id', 'element_type', 'nnodes'])
            size_end = n

            # max sizes
            self._size_node_start = 0
            self._size_node_end = n
            self._size_element_start = 0
            self._size_element_end = nelements_size

            self._inode_start = 0
            self._inode_end = nnodes
            self._ielement_start = 0
            self._ielement_end = nelements
            #print('_inode=%s _ielement=%s' %(nnodes, nelements))
        return (self._inode_start, self._inode_end, self._ielement_start, self._ielement_end)

    def _finalize(self, dt):
        ndt, nnodes, nelements, dts = self._get_shape()

        if dt is not None and dt != dts[-1]:
            return
        #print("----finalize----")

        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append('C' if grid_type==0 else grid_type)
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        update_index = True
        #print("final A\n", self.data.to_string())
        if update_index:
            if dts[0] is not None:
                name = self.data_code['name']
                self.data = self.data.set_index([name, 'element_id', 'node_id'])
            else:
                self.data = self.data.set_index(['element_id', 'node_id'])
            #self.element_data = self.element_data.set_index(['element_id'])

        #print(self.data.to_string())
        #print("final\n", self.data.to_string())
        #print(self.element_data.to_string())
        del self._inode_start
        del self._inode_end
        #print('---PlateStressObject---')
        #print(self.data.to_string())

    def _is_full(self, nnodes, nelements):
        self._inode_start += nnodes
        self._ielement_start += nelements
        
        #print('check....inode_start=%r size_node_end=%r' %(self._inode_start, self._size_node_end))
        #print('check...._ielement_start=%s size_element_end=%s' %(self._ielement_start, self._size_element_end))
        if self._inode_start == self._size_node_end:
            return True
        elif self._ielement_start == self._size_element_end:
            self._ielement_start = 0
        return False

    def get_stats(self):
        msg = self._get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.radial)
            r0 = self.radial.keys()[0]
            nelements = len(self.radial[r0])
            dt_string = name + ', '
            msg.append('  real type=%s n%ss=%s nelements=%s\n'
                       % (self.__class__.__name__, name, ntimes, nelements))
        else:
            nelements = len(self.radial)
            dt_string = ''
            msg.append('  real type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        #msg.append('  eType, radial, azimuthal, axial, shear, '
                   #'omax, oms, ovm\n')
        msg.append('  data: index  : %selement_id, node_id\n' % dt_string)
        msg.append('      : results: radial, azimuthal, axial, shear, omax, oms, ovm\n')
        return msg


class TriaxStressObject(RealTriaxObject, StressObject):
    """
    ::

      # format_code=1 sort_code=0 stressCode=0
                                        S T R E S S E S   I N   T R I A X 6   E L E M E N T S
      ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES
         ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR
         5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                  4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
    """
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        RealTriaxObject.__init__(self)
        StressObject.__init__(self, data_code, isubcase, read_mode)
        self.eType = 'CTRIAX6'

        #self.code = [self.format_code, self.sort_code, self.s_code]
        self.radial = {}
        self.azimuthal = {}
        self.axial = {}
        self.shear = {}
        self.omax = {}
        self.oms = {}
        self.ovm = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def _get_headers(self):
        return ['radial', 'azimuthal', 'axial', 'shear', 'omax', 'oms', 'ovm']

    def add_f06_data(self, data, transient):
        raise Exception('Not Implemented')
        if transient is None:
            for line in data:
                (eid, axial, MSa, torsion, MSt) = line
                if MSa is None:
                    MSa = 0.
                if MSt is None:
                    MSt = 0.
                self.axial[eid] = axial
                self.MS_axial[eid] = MSa
                self.torsion[eid] = torsion
                self.MS_torsion[eid] = MSt
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.s1:
            self.update_dt(self.data_code, dt)
            self.isTransient = True

        for line in data:
            (eid, axial, MSa, torsion, MSt) = line
            if MSa is None:
                MSa = 0.
            if MSt is None:
                MSt = 0.
            self.axial[dt][eid] = axial
            self.MS_axial[dt][eid] = MSa
            self.torsion[dt][eid] = torsion
            self.MS_torsion[dt][eid] = MSt

    def get_transients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.radial[dt] = {}
        self.azimuthal[dt] = {}
        self.axial[dt] = {}
        self.shear[dt] = {}
        self.omax[dt] = {}
        self.oms[dt] = {}
        self.ovm[dt] = {}

    def add_new_eid(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #print "**?eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid] = {nid: rs}
        self.azimuthal[eid] = {nid: azs}
        self.axial[eid] = {nid: As}
        self.shear[eid] = {nid: ss}
        self.omax[eid] = {nid: maxp}
        self.oms[eid] = {nid: tmax}
        self.ovm[eid] = {nid: octs}

    def add(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid][nid] = rs
        self.azimuthal[eid][nid] = azs
        self.axial[eid][nid] = As
        self.shear[eid][nid] = ss
        self.omax[eid][nid] = maxp
        self.oms[eid][nid] = tmax
        self.ovm[eid][nid] = octs

    def add_new_eid_sort1(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #assert isinstance(eid,int)
        #assert eid >= 0
        #print "*  eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        if dt not in self.radial:
            self.add_new_transient(dt)
        self.radial[dt][eid] = {nid: rs}
        self.azimuthal[dt][eid] = {nid: azs}
        self.axial[dt][eid] = {nid: As}
        self.shear[dt][eid] = {nid: ss}
        self.omax[dt][eid] = {nid: maxp}
        self.oms[dt][eid] = {nid: tmax}
        self.ovm[dt][eid] = {nid: octs}

    def add_sort1(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[dt][eid][nid] = rs
        self.azimuthal[dt][eid][nid] = azs
        self.axial[dt][eid][nid] = As
        self.shear[dt][eid][nid] = ss
        self.omax[dt][eid][nid] = maxp
        self.oms[dt][eid][nid] = tmax
        self.ovm[dt][eid][nid] = octs

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        msg = header + ['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
                        '   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                        '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        #out = []
        for eid, radial in sorted(self.radial.iteritems()):
            for nid in sorted(radial):
                rad = self.radial[eid][nid]
                azimuth = self.azimuthal[eid][nid]
                axial = self.axial[eid][nid]
                shear = self.shear[eid][nid]
                omax = self.omax[eid][nid]
                oms = self.oms[eid][nid]
                ovm = self.ovm[eid][nid]
                if nid == 0:
                    Eid = eid
                else:
                    Eid = ''
                ([rad, azimuth, axial, shear, omax, oms, ovm], isAllZeros) = writeFloats13E([rad, azimuth, axial, shear, omax, oms, ovm])
                msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, radial, azimuth, axial, shear, omax, oms, ovm.rstrip()))
            msg.append('\n')

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, f,
                          pageNum=1, is_mag_phase=False):
        words = ['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
                 '   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                 '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        msg = []
        for dt, Radial in sorted(self.radial.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, radial in sorted(Radial.iteritems()):
                for nid in sorted(radial):
                    rad = self.radial[dt][eid][nid]
                    azimuth = self.azimuthal[dt][eid][nid]
                    axial = self.axial[dt][eid][nid]
                    shear = self.shear[dt][eid][nid]
                    omax = self.omax[dt][eid][nid]
                    oms = self.oms[dt][eid][nid]
                    ovm = self.ovm[dt][eid][nid]
                    if nid == 0:
                        Eid = eid
                    else:
                        Eid = ''
                    ([rad, azimuth, axial, shear, omax, oms, ovm], isAllZeros) = writeFloats13E([rad, azimuth, axial, shear, omax, oms, ovm])
                    msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, rad, azimuth, axial, shear, omax, oms, ovm.rstrip()))
                msg.append('\n')

            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1

class TriaxStrainObject(RealTriaxObject, StrainObject):
    """
    ::

      # format_code=1 sort_code=0 stressCode=0
                                        S T R A I N S   I N   T R I A X 6   E L E M E N T S
      ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES
         ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR
         5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                  4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
    """
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        RealTriaxObject.__init__(self)
        StrainObject.__init__(self, data_code, isubcase, read_mode)
        self.eType = 'CTRIAX6'

        #self.code = [self.format_code, self.sort_code, self.s_code]
        self.radial = {}
        self.azimuthal = {}
        self.axial = {}
        self.shear = {}
        self.emax = {}
        self.ems = {}
        self.evm = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def _get_headers(self):
        return ['radial', 'azimuthal', 'axial', 'shear', 'emax', 'ems', 'evm']

    def add_f06_data(self, data, transient):
        raise Exception('Not Implemented')

    def get_transients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.radial[dt] = {}
        self.azimuthal[dt] = {}
        self.axial[dt] = {}
        self.shear[dt] = {}
        self.emax[dt] = {}
        self.ems[dt] = {}
        self.evm[dt] = {}

    def add_new_eid(self, dt, eid, nid, rs, azs, As, ss, maxp, tmax, octs):
        #print "**?eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid] = {nid: rs}
        self.azimuthal[eid] = {nid: azs}
        self.axial[eid] = {nid: As}
        self.shear[eid] = {nid: ss}
        self.emax[eid] = {nid: maxp}
        self.ems[eid] = {nid: emax}
        self.evm[eid] = {nid: ects}

    def add(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid][nid] = rs
        self.azimuthal[eid][nid] = azs
        self.axial[eid][nid] = As
        self.shear[eid][nid] = ss
        self.emax[eid][nid] = maxp
        self.ems[eid][nid] = emax
        self.evm[eid][nid] = ects

    def add_new_eid_sort1(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        #assert isinstance(eid,int)
        #assert eid >= 0
        #print "*  eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[dt][eid] = {nid: rs}
        self.azimuthal[dt][eid] = {nid: azs}
        self.axial[dt][eid] = {nid: As}
        self.shear[dt][eid] = {nid: ss}
        self.emax[dt][eid] = {nid: maxp}
        self.ems[dt][eid] = {nid: emax}
        self.evm[dt][eid] = {nid: ects}

    def add_sort1(self, dt, eid, nid, rs, azs, As, ss, maxp, emax, ects):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[dt][eid][nid] = rs
        self.azimuthal[dt][eid][nid] = azs
        self.axial[dt][eid][nid] = As
        self.shear[dt][eid][nid] = ss
        self.emax[dt][eid][nid] = maxp
        self.ems[dt][eid][nid] = emax
        self.evm[dt][eid][nid] = ects

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        msg = header + ['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
                        '   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                        '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        #out = []
        for eid, radial in sorted(self.radial.iteritems()):
            for nid in sorted(radial):
                rad = self.radial[eid][nid]
                azimuth = self.azimuthal[eid][nid]
                axial = self.axial[eid][nid]
                shear = self.shear[eid][nid]
                emax = self.emax[eid][nid]
                ems = self.ems[eid][nid]
                evm = self.evm[eid][nid]
                if nid == 0:
                    Eid = eid
                else:
                    Eid = ''
                ([rad, azimuth, axial, shear, emax, ems, evm], isAllZeros) = writeFloats13E([rad, azimuth, axial, shear, emax, ems, evm])
                msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' % (Eid, nid, radial, azimuth, axial, shear, emax, ems, evm.rstrip()))
            msg.append('\n')

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, f,
                          pageNum=1, is_mag_phase=False):
        words = ['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
                 '   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
                 '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n', ]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        msg = []
        for dt, Radial in sorted(self.radial.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, radial in sorted(Radial.iteritems()):
                for nid in sorted(radial):
                    rad = self.radial[dt][eid][nid]
                    azimuth = self.azimuthal[dt][eid][nid]
                    axial = self.axial[dt][eid][nid]
                    shear = self.shear[dt][eid][nid]
                    emax = self.emax[dt][eid][nid]
                    ems = self.ems[dt][eid][nid]
                    evm = self.evm[dt][eid][nid]
                    if nid == 0:
                        Eid = eid
                    else:
                        Eid = ''
                    ([rad, azimuth, axial, shear, emax, ems, evm], isAllZeros) = writeFloats13E([rad, azimuth, axial, shear, emax, ems, evm])
                    msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n'
                               % (Eid, nid, rad, azimuth, axial, shear, emax,
                                  ems, evm.rstrip()))
                msg.append('\n')

            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1