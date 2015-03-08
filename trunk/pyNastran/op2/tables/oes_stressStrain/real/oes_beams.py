from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from itertools import count

from numpy import zeros
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E


class RealBeamArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        if is_sort1:
            pass
            #if dt is not None:
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
                #self.addNewNode = self.addNewNodeSort1
        else:
            raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2
            #self.addNewNode = self.addNewNodeSort2

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
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 2:
            nnodes_per_element = 10
        else:
            raise NotImplementedError(self.element_type)

        self.nnodes = nnodes_per_element
        self.nelements //= self.ntimes
        self.ntotal = self.nelements  #* 2  # for A/B
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

        # sxc, sxd, sxe, sxf
        # smax, smin, MSt, MSc
        self.xxb = zeros(self.ntotal, dtype='float32')
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def add_new_eid(self, dt, eid, out):
        self.add_new_eid_sort1(dt, eid, out)

    def add_new_eid_sort1(self, dt, eid, out):
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        assert isinstance(eid, int), eid
        assert eid >= 0, eid
        self._times[self.itime] = dt
        self.element_node[self.itotal] = [eid, grid]
        self.xxb[self.itotal] = sd
        self.data[self.itime, self.itotal, :] = [sxc, sxd, sxe, sxf,
                                                 smax, smin, mst, msc]
        self.itotal += 1
        self.ielement += 1

    def add(self, dt, eid, out):
        self.add_sort1(dt, eid, out)

    def add_sort1(self, dt, eid, out):
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out

        self.element_node[self.itotal, :] = [eid, grid]
        self.xxb[self.itotal] = sd
        self.data[self.itime, self.itotal, :] = [sxc, sxd, sxe, sxf,
                                                 smax, smin, mst, msc]
        self.itotal += 1

    #def add_sort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert eid is not None
        #msg = "i=%s dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (
            #self.itotal, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        ##print(msg)
        #if isinstance(nodeID, string_types):
            #nodeID = 0
        ##assert isinstance(nodeID, int), nodeID
        #self.element_node[self.itotal, :] = [eid, nodeID]
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy, txy, angle, majorP, minorP, ovm]
        #self.itotal += 1

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
        #nlayers = 2
        nelements = self.ntotal // self.nnodes  # // 2

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes_per_element=%i ntotal=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, ntotal))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i nnodes_per_element=%i ntotal=%i\n'
                       % (self.__class__.__name__, nelements, nnodes, ntotal))
            ntimes_word = 1
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element types: %s\n  ' % ', '.join(self.element_names))
        msg += self.get_data_code()
        return msg

    #def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        #return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        msg = self._get_msgs()
        (ntimes, ntotal, four) = self.data.shape

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        xxbs = self.xxb
        #print('CBEAM ntimes=%s ntotal=%s' % (ntimes, ntotal))
        for itime in range(ntimes):
            dt = self._times[itime]
            if self.nonlinear_factor is not None:
                dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                header[1] = dtLine
                if hasattr(self, 'eigr'):
                    header[2] = ' %14s = %12.6E\n' % ('EIGENVALUE', self.eigrs[itime])
            f.write(''.join(header + msg))

            sxcs = self.data[itime, :, 0]
            sxds = self.data[itime, :, 1]
            sxes = self.data[itime, :, 2]
            sxfs = self.data[itime, :, 3]
            smaxs = self.data[itime, :, 4]
            smins = self.data[itime, :, 5]
            SMts = self.data[itime, :, 6]
            SMcs = self.data[itime, :, 7]

            # loop over all the elements
            eid_old = None
            xxb_old = None
            for (i, eid, nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) in zip(
                count(), eids, nids, xxbs, sxcs, sxds, sxes, sxfs, smaxs, smins, SMts, SMcs):
                if eid != eid_old:
                    f.write('0  %8i\n' % eid)
                if xxb == xxb_old:
                    continue
                # #if eid != eid_old and xxb != xxb_old:
                    #continue
                #asdf
                vals = [sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc] = vals2
                f.write('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf,
                                                                                    sMax, sMin, SMt, SMc.strip()))
                eid_old = eid
                xxb_old = xxb

            f.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num


class RealBeamStressArray(RealBeamArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBeamArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def isStress(self):
        return True
    def isStrain(self):
        return False

    def get_headers(self):
        headers = [
            #'grid', 'xxb',
            'sxc', 'sxd', 'sxe', 'sxf',
            'smax', 'smin', 'MS_tension', 'MS_compression'
        ]
        return headers

    def _get_msgs(self):
        if self.element_type == 2:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = ['                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/\n',
                        '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        return msg


class RealBeamStrainArray(RealBeamArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBeamArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def isStress(self):
        return False
    def isStrain(self):
        return True

    def get_headers(self):
        headers = [
            #'grid', 'xxb',
            'sxc', 'sxd', 'sxe', 'sxf',
            'smax', 'smin', 'MS_tension', 'MS_compression'
        ]
        return headers

    def _get_msgs(self):
        if self.element_type == 2:
            pass
        else:
            raise NotImplementedError(self.element_type)
        msg = ['                                  S T R A I N S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/\n',
                        '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        return msg


class RealBeamStress(StressObject):
    """
    ::

      [1,0,0]
                   S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )
                        STAT DIST/
       ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C
              1       1   0.000   -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04
                      2   1.000   -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CBEAM'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.xxb = {}
        self.grids = {}
        self.smax = {}
        self.smin = {}
        self.MS_tension = {}
        self.MS_compression = {}

        self.sxc = {}
        self.sxd = {}
        self.sxe = {}
        self.sxf = {}
        #self.isImaginary = False

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        #else:
            #assert dt is not None
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.smax)
            s0 = self.smax.keys()[0]
            nelements = len(self.smax[s0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.smax)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, xxb, grids, smax, smin, MS_tension, '
                   'MS_compression, sxc, sxd, sxe, sxf\n')
        return msg

    def getLengthTotal(self):
        return 444  # 44+10*40   (11 nodes)

    def getLength1(self):
        return (44, 'ifffffffff')

    def getLength2(self):
        return (40, 'ifffffffff')

    def delete_transient(self, dt):
        del self.sxc[dt]
        del self.sxd[dt]
        del self.sxe[dt]
        del self.sxf[dt]
        del self.smax[dt]
        del self.smin[dt]
        del self.MS_tension[dt]
        del self.MS_compression[dt]

    def get_transients(self):
        k = self.smax.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #print("addNewTransient_beam+1+0")
        self.dt = dt
        self.sxc[dt] = {}
        self.sxd[dt] = {}
        self.sxe[dt] = {}
        self.sxf[dt] = {}
        self.smax[dt] = {}
        self.smin[dt] = {}
        self.MS_tension[dt] = {}
        self.MS_compression[dt] = {}

    def add_new_eid(self, dt, eid, out):
        #print("Beam Stress add_new_eid...")
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        #print("eid=%s grid=%s" % (eid, grid))
        assert eid >= 0
        #assert isinstance(eid, int)
        #assert isinstance(grid, int)
        self.grids[eid] = [grid]
        self.xxb[eid] = [sd]
        self.sxc[eid] = [sxc]
        self.sxd[eid] = [sxd]
        self.sxe[eid] = [sxe]
        self.sxf[eid] = [sxf]
        self.smax[eid] = [smax]
        self.smin[eid] = [smin]
        self.MS_tension[eid] = [mst]
        self.MS_compression[eid] = [msc]
        return eid

    def add_f06_data(self, data, dt):
        if dt is not None:
            raise NotImplementedError(dt)
        for datai in data:
            (eid, grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = datai
            if eid in self.grids:
                self.grids[eid].append(grid)
                self.xxb[eid].append(sd)
                self.sxc[eid].append(sxc)
                self.sxd[eid].append(sxd)
                self.sxe[eid].append(sxe)
                self.sxf[eid].append(sxf)
                self.smax[eid].append(smax)
                self.smin[eid].append(smin)
                self.MS_tension[eid].append(mst)
                self.MS_compression[eid].append(msc)
            else:
                self.grids[eid] = [grid]
                self.xxb[eid] = [sd]
                self.sxc[eid] = [sxc]
                self.sxd[eid] = [sxd]
                self.sxe[eid] = [sxe]
                self.sxf[eid] = [sxf]
                self.smax[eid] = [smax]
                self.smin[eid] = [smin]
                self.MS_tension[eid] = [mst]
                self.MS_compression[eid] = [msc]

    def add_new_eid_sort1(self, dt, eid, out):
        #print "Beam Transient Stress add_new_eid..."
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out

        assert eid >= 0
        if dt not in self.sxc:
            self.add_new_transient(dt)
        self.grids[eid] = [grid]
        self.xxb[eid] = [sd]
        self.sxc[dt][eid] = [sxc]
        self.sxd[dt][eid] = [sxd]
        self.sxe[dt][eid] = [sxe]
        self.sxf[dt][eid] = [sxf]
        self.smax[dt][eid] = [smax]
        self.smin[dt][eid] = [smin]
        self.MS_tension[dt][eid] = [mst]
        self.MS_compression[dt][eid] = [msc]
        return eid

    def add(self, dt, eid, out):
        #print "Beam Stress add..."
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        if grid:
            self.grids[eid].append(grid)
            self.xxb[eid].append(sd)
            self.sxc[eid].append(sxc)
            self.sxd[eid].append(sxd)
            self.sxe[eid].append(sxe)
            self.sxf[eid].append(sxf)
            self.smax[eid].append(smax)
            self.smin[eid].append(smin)
            self.MS_tension[eid].append(mst)
            self.MS_compression[eid].append(msc)

    def add_sort1(self, dt, eid, out):
        #print "Beam Transient Stress add..."
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        if grid:
            self.grids[eid].append(grid)
            self.xxb[eid].append(sd)
            #self.sd[dt][eid].append(sd)
            self.sxc[dt][eid].append(sxc)
            self.sxd[dt][eid].append(sxd)
            self.sxe[dt][eid].append(sxe)
            self.sxf[dt][eid].append(sxf)
            self.smax[dt][eid].append(smax)
            self.smin[dt][eid].append(smin)
            self.MS_tension[dt][eid].append(mst)
            self.MS_compression[dt][eid].append(msc)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f)

        msg = header + ['                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/\n',
                        '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']

        for eid in sorted(self.smax):
            msg.append('0  %8i\n' % (eid))
            #print self.xxb[eid]
            for i, nid in enumerate(self.grids[eid]):
                #print i,nid
                xxb = self.xxb[eid][i]
                sxc = self.sxc[eid][i]
                sxd = self.sxd[eid][i]
                sxe = self.sxe[eid][i]
                sxf = self.sxf[eid][i]
                sMax = self.smax[eid][i]
                sMin = self.smin[eid][i]
                SMt = self.MS_tension[eid][i]
                SMc = self.MS_compression[eid][i]
                (vals2, is_all_zeros) = writeFloats13E([sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc])
                (sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) = vals2
                msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc.strip()))

        msg.append(pageStamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/\n',
                 '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        msg = []
        for dt, SMaxs in sorted(iteritems(self.smax)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, Smax in sorted(iteritems(SMaxs)):
                msg.append('0  %8i\n' % (eid))
                for i, nid in enumerate(self.grids[eid]):
                    xxb = self.xxb[eid][i]
                    #sd  = self.sd[eid][i]
                    sxc = self.sxc[dt][eid][i]
                    sxd = self.sxd[dt][eid][i]
                    sxe = self.sxe[dt][eid][i]
                    sxf = self.sxf[dt][eid][i]
                    sMax = self.smax[dt][eid][i]
                    sMin = self.smin[dt][eid][i]
                    SMt = self.MS_tension[dt][eid][i]
                    SMc = self.MS_compression[dt][eid][i]
                    (vals2, is_all_zeros) = writeFloats13E([sxc, sxd,
                                                          sxe, sxf, sMax, sMin, SMt, SMc])
                    (sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) = vals2
                    msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc.strip()))

            msg.append(pageStamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1


class RealBeamStrain(StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = 'CBEAM'  # {} # 'CBEAM/CONBEAM'

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.xxb = {}
        self.grids = {}
        self.sxc = {}
        self.sxd = {}
        self.sxe = {}
        self.sxf = {}
        self.smax = {}
        self.smin = {}
        self.MS_tension = {}
        self.MS_compression = {}

        if is_sort1:
            if dt is not None:
                #self.add_new_transient = self.add_new_transient
                self.add_new_eid = self.add_new_eid_sort1
                self.add = self.add_sort1
                self.add_new_transient(dt)
            #self.add_new_eid = self.add_new_eid
            #self.add = self.add
        #else:
            #raise  RuntimeError("Invalid Code: beamStress - get the format/sort/stressCode=%s" % (self.code))
        #if dt is not None:
            #self.isTransient = True
            #self.dt = self.nonlinear_factor
            #self.add_new_transient(dt)

    def add_f06_data(self, data, dt):
        if dt is not None:
            raise NotImplementedError(dt)
        for datai in data:
            #print('smax', self.smax)
            (eid, grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = datai
            #print('(*')
            if eid in self.grids:
                self.grids[eid].append(grid)
                self.xxb[eid].append(sd)
                self.sxc[eid].append(sxc)
                self.sxd[eid].append(sxd)
                self.sxe[eid].append(sxe)
                self.sxf[eid].append(sxf)
                self.smax[eid].append(smax)
                self.smin[eid].append(smin)
                self.MS_tension[eid].append(mst)
                self.MS_compression[eid].append(msc)
            else:
                self.grids[eid] = [grid]
                self.xxb[eid] = [sd]
                self.sxc[eid] = [sxc]
                self.sxd[eid] = [sxd]
                self.sxe[eid] = [sxe]
                self.sxf[eid] = [sxf]
                self.smax[eid] = [smax]
                self.smin[eid] = [smin]
                self.MS_tension[eid] = [mst]
                self.MS_compression[eid] = [msc]
        #print('smax', self.smax)

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.smax)
            s0 = self.smax.keys()[0]
            nelements = len(self.smax[s0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.smax)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, xxb, grids, smax, smin, MS_tension, '
                   'MS_compression, sxc, sxd, sxe, sxf\n')
        return msg

    def getLengthTotal(self):
        return 444  # 44+10*40   (11 nodes)

    def getLength1(self):
        return (44, 'i9f')

    def getLength2(self):
        return (40, 'i9f')

    def delete_transient(self, dt):
        del self.sxc[dt]
        del self.sxd[dt]
        del self.sxe[dt]
        del self.sxf[dt]
        del self.smax[dt]
        del self.smin[dt]
        del self.MS_tension[dt]
        del self.MS_compression[dt]

    def get_transients(self):
        k = self.smax.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        .. note:: make sure you set self.dt first
        """
        #print "addNewTransient_beam+1+0"
        self.dt = dt
        self.grids = {}
        self.xxb = {}
        self.sxc[dt] = {}
        self.sxd[dt] = {}
        self.sxe[dt] = {}
        self.sxf[dt] = {}
        self.smax[dt] = {}
        self.smin[dt] = {}
        self.MS_tension[dt] = {}
        self.MS_compression[dt] = {}

    def add_new_eid(self, dt, eid, out):
        #print "Beam Stress add_new_eid..."
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        #print "eid=%s grid=%s" %(eid,grid)
        assert eid >= 0
        #assert isinstance(eid,int)
        #assert isinstance(grid,int)
        self.grids[eid] = [grid]
        self.xxb[eid] = [sd]
        self.sxc[eid] = [sxc]
        self.sxd[eid] = [sxd]
        self.sxe[eid] = [sxe]
        self.sxf[eid] = [sxf]
        self.smax[eid] = [smax]
        self.smin[eid] = [smin]
        self.MS_tension[eid] = [mst]
        self.MS_compression[eid] = [msc]
        return eid

    def add(self, dt, eid, out):
        #print "Beam Stress add..."
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        if grid:
            self.grids[eid].append(grid)
            self.xxb[eid].append(sd)
            self.sxc[eid].append(sxc)
            self.sxd[eid].append(sxd)
            self.sxe[eid].append(sxe)
            self.sxf[eid].append(sxf)
            self.smax[eid].append(smax)
            self.smin[eid].append(smin)
            self.MS_tension[eid].append(mst)
            self.MS_compression[eid].append(msc)

    def add_new_eid_sort1(self, dt, eid, out):
        #print "Beam Transient Stress add_new_eid..."
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out

        assert eid >= 0
        if dt not in self.sxc:
            self.add_new_transient(dt)

        self.grids[eid] = [grid]
        self.xxb[eid] = [sd]
        self.sxc[dt][eid] = [sxc]
        self.sxd[dt][eid] = [sxd]
        self.sxe[dt][eid] = [sxe]
        self.sxf[dt][eid] = [sxf]
        self.smax[dt][eid] = [smax]
        self.smin[dt][eid] = [smin]
        self.MS_tension[dt][eid] = [mst]
        self.MS_compression[dt][eid] = [msc]
        return eid

    def add_sort1(self, dt, eid, out):
        (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        if grid:
            self.grids[eid].append(grid)
            self.xxb[eid].append(sd)
            self.sxc[dt][eid].append(sxc)
            self.sxd[dt][eid].append(sxd)
            self.sxe[dt][eid].append(sxe)
            self.sxf[dt][eid].append(sxf)
            self.smax[dt][eid].append(smax)
            self.smin[dt][eid].append(smin)
            self.MS_tension[dt][eid].append(mst)
            self.MS_compression[dt][eid].append(msc)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f)

        msg = header + ['                                  S T R A I N S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/\n',
                        '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']

        for eid in sorted(self.smax):
            msg.append('0  %8i\n' % eid)
            #print self.xxb[eid]
            #print("self.grids =", self.grids)
            #print("eid =", eid)
            for i, nid in enumerate(self.grids[eid]):
                xxb = self.xxb[eid][i]
                sxc = self.sxc[eid][i]
                sxd = self.sxd[eid][i]
                sxe = self.sxe[eid][i]
                sxf = self.sxf[eid][i]
                sMax = self.smax[eid][i]
                sMin = self.smin[eid][i]
                SMt = self.MS_tension[eid][i]
                SMc = self.MS_compression[eid][i]
                (vals2, is_all_zeros) = writeFloats13E([
                    sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc])
                (sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) = vals2
                msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc))

        msg.append(pageStamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                  S T R A I N S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/\n',
                 '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        msg = []
        for dt, SMaxs in sorted(iteritems(self.smax)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, Smax in sorted(iteritems(SMaxs)):
                msg.append('0  %8i\n' % (eid))
                for i, nid in enumerate(self.grids[eid]):
                    xxb = self.xxb[eid][i]
                    sxc = self.sxc[dt][eid][i]
                    sxd = self.sxd[dt][eid][i]
                    sxe = self.sxe[dt][eid][i]
                    sxf = self.sxf[dt][eid][i]
                    sMax = self.smax[dt][eid][i]
                    sMin = self.smin[dt][eid][i]
                    SMt = self.MS_tension[dt][eid][i]
                    SMc = self.MS_compression[dt][eid][i]
                    (vals2, is_all_zeros) = writeFloats13E([sxc, sxd, sxe, sxf,
                                                          sMax, sMin, SMt, SMc])
                    (sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) = vals2
                    msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc))

            msg.append(pageStamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1
