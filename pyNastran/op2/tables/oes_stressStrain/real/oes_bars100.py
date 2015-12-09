from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems

from itertools import count
from numpy import zeros, searchsorted, ravel
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header


class RealBar10NodesArray(OES_Object):
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
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                #self.addNewNode = self.addNewNodeSort1
        else:
            raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.add_sort2
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
        # print('RealBar10NodesArray isubcase=%s ntimes=%s nelements=%s ntotal=%s' % (
            # self.isubcase, self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type == 100:
            nnodes_per_element = 1
        else:
            raise NotImplementedError(self.element_type)

        self.nnodes = nnodes_per_element
        self.nelements //= self.ntimes
        #self.ntotal = self.nelements  #* 2  # for A/B
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
        self.element = zeros(self.ntotal, dtype='int32')

        #[sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]
        self.data = zeros((self.ntimes, self.ntotal, 9), dtype='float32')

    def add_new_eid(self, eType, dt, eid, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS):
        self.add_new_eid_sort1(eType, dt, eid,
                               sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS)

    def add_new_eid_sort1(self, eType, dt, eid,
                          sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS):
        self._times[self.itime] = dt
        # print('isubcase=%s itotal=%s ieid=%s eid=%s' % (self.isubcase, self.itotal, self.ielement, eid))
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]
        self.itotal += 1
        self.ielement += 1

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
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ##ind = searchsorted(eids, self.element)
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = self._get_msgs()
        #print('CBAR ntimes=%s ntotal=%s' % (ntimes, ntotal))
        if self.is_sort1():
            page_num = self._write_sort1_as_sort1(f, header, page_stamp, msg, page_num)
        else:
            raise RuntimeError()
        return page_num

    def _write_sort1_as_sort1(self, f06_file, header, page_stamp, msg, page_num):
        (ntimes, ntotal) = self.data.shape[:2]
        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            sd = self.data[itime, :, 0]
            sxc = self.data[itime, :, 1]
            sxd = self.data[itime, :, 2]
            sxe = self.data[itime, :, 3]
            sxf = self.data[itime, :, 4]
            axial = self.data[itime, :, 5]
            smax = self.data[itime, :, 6]
            smin = self.data[itime, :, 7]
            MS = self.data[itime, :, 8]

            for (i, eid, sdi, sxci, sxdi, sxei, sxfi, axiali, smaxi, smini, MSi) in zip(
                     count(), eids, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS):

                vals = [sdi, sxci, sxdi, sxei, sxfi, axiali, smaxi, smini, MSi]
                vals2 = write_floats_13e(vals)
                [sdi, sxci, sxdi, sxei, sxfi, axiali, smaxi, smini, MSi] = vals2
                f06_file.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s %s\n'
                               % (eid, sdi, sxci, sxdi, sxei, sxfi, axiali, smaxi, smini, MSi))

            f06_file.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num

class RealBar10NodesStressArray(RealBar10NodesArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBar10NodesArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        #if self.is_fiber_distance():
            #fiber_dist = 'fiber_distance'
        #else:
            #fiber_dist = 'fiber_curvature'

        #if self.is_von_mises():
            #ovm = 'von_mises'
        #else:
            #ovm = 'max_shear'
        headers = ['sd', 'sxc', 'sxd', 'sxe', 'sxf', 'axial', 'smax', 'smin', 'MS']
        return headers

    def _get_msgs(self):
        msg = [
            '                         S T R E S S   D I S T R I B U T I O N   I N   B A R   E L E M E N T S       ( C B A R )\n'
            '0    ELEMENT  STATION    SXC           SXD           SXE           SXF            AXIAL          S-MAX         S-MIN         M.S.-T\n'
            '       ID.     (PCT)                                                                                                         M.S.-C\n'
            #'            1   0.000   4.919032E+05 -4.348710E+05 -4.348710E+05  4.919032E+05   0.0            4.919032E+05 -4.348710E+05 \n'
        ]
        return msg


class RealBar10NodesStrainArray(RealBar10NodesArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBar10NodesArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        #if self.is_fiber_distance():
            #fiber_dist = 'fiber_distance'
        #else:
            #fiber_dist = 'fiber_curvature'

        #if self.is_von_mises():
            #ovm = 'von_mises'
        #else:
            #ovm = 'max_shear'
        headers = ['sd', 'sxc', 'sxd', 'sxe', 'sxf', 'axial', 'smax', 'smin', 'MS']
        return headers

    def _get_msgs(self):
        msg = [
            '                         S T R A I N   D I S T R I B U T I O N   I N   B A R   E L E M E N T S       ( C B A R )\n'
            '0    ELEMENT  STATION    SXC           SXD           SXE           SXF            AXIAL          S-MAX         S-MIN         M.S.-T\n'
            '       ID.     (PCT)                                                                                                         M.S.-C\n'
            #'            1   0.000   4.919032E+05 -4.348710E+05 -4.348710E+05  4.919032E+05   0.0            4.919032E+05 -4.348710E+05 \n'
        ]
        return msg

