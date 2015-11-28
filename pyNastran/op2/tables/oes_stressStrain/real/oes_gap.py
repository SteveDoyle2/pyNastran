from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import count
from numpy import zeros, searchsorted, ravel

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, _eigenvalue_header


class NonlinearGapStressArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        #print(self.code_information())

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
        headers = ['compX', 'shearY', 'shearZ', 'axialU', 'shearV', 'shearW', 'slipV', 'slipW']
        return headers

    def build(self):
        if self.is_built:
            return
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        #if self.element_type == None?:
            #nnodes_per_element = 1
        #else:
            #raise NotImplementedError(self.element_type)

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

        # [compX, shearY, shearZ, axialU, shearV, shearW, slipV, slipW]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def add_sort1(self, dt, eid, compX, shearY, shearZ, axialU, shearV, shearW, slipV, slipW, form1, form2):
        assert isinstance(eid, int)
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [compX, shearY, shearZ, axialU, shearV, shearW, slipV, slipW]
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
        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            compX, shearY, shearZ, axialU, shearV, shearW, slipV, slipW
            compX = self.data[itime, :, 0]
            shearY = self.data[itime, :, 1]
            shearZ = self.data[itime, :, 2]
            axialU = self.data[itime, :, 3]
            shearV = self.data[itime, :, 4]
            shearW = self.data[itime, :, 5]
            slipV = self.data[itime, :, 6]
            slipW = self.data[itime, :, 7]

            # loop over all the elements
            for (i, eid, compXi, shearYi, shearZi, axialUi, shearVi, shearWi, slipVi, slipWi) in zip(
                count(), eids, compX, shearY, shearZ, axialU, shearV, shearW, slipV, slipW):

                vals = [compXi, shearYi, shearZi, axialUi, shearVi, shearWi, slipVi, slipWi]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [compXi, shearYi, shearZi, axialUi, shearVi, shearWi, slipVi, slipWi] = vals2
                f.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                    % (eid, compXi, shearYi, shearZi, axialUi, shearVi, shearWi, slipVi, slipWi))
            f.write(page_stamp % page_num)
            page_num += 1
        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num


class NonlinearGapStress(StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.element_name = self.data_code['element_name']

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.compX = {}
        self.shearY = {}
        self.shearZ = {}
        self.axialU = {}
        self.shearV = {}
        self.shearW = {}
        self.slipV = {}
        self.slipW = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.add_sort2
            self.add_new_eid = self.add_new_eid_sort2

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def get_stats(self):
        nelements = len(self.eType)
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.compX)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, compX, shearY, shearZ, axialU, shearV, shearW, slipV, slipW\n')
        msg.append('  %s\n' % self.element_name)
        return msg

    def getLength(self):
        return (8, 'f')

    def delete_transient(self, dt):
        del self.compX[dt]

    def get_transients(self):
        k = self.compX.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.element_name = self.data_code['element_name']
        self.dt = dt
        self.compX[dt] = {}
        self.shearY[dt] = {}
        self.shearZ[dt] = {}
        self.axialU[dt] = {}
        self.shearV[dt] = {}
        self.shearW[dt] = {}
        self.slipV[dt] = {}
        self.slipW[dt] = {}

    def add_new_eid(self, dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        self.eType[eid] = self.element_name
        self.compX[eid] = cpx
        self.shearY[eid] = shy
        self.shearZ[eid] = shz
        self.axialU[eid] = au
        self.shearV[eid] = shv
        self.shearW[eid] = shw
        self.slipV[eid] = slv
        self.slipW[eid] = slp

    def add_new_eid_sort1(self, dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        if dt not in self.compX:
            self.add_new_transient(dt)
        self.eType[eid] = self.element_name
        self.compX[dt][eid] = cpx
        self.shearY[dt][eid] = shy
        self.shearZ[dt][eid] = shz
        self.axialU[dt][eid] = au
        self.shearV[dt][eid] = shv
        self.shearW[dt][eid] = shw
        self.slipV[dt][eid] = slv
        self.slipW[dt][eid] = slp

    def add_new_eid_sort2(self, eid, dt, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        if dt not in self.compX:
            self.add_new_transient(dt)
        self.eType[eid] = self.element_name
        self.compX[dt][eid] = cpx
        self.shearY[dt][eid] = shy
        self.shearZ[dt][eid] = shz
        self.axialU[dt][eid] = au
        self.shearV[dt][eid] = shv
        self.shearW[dt][eid] = shw
        self.slipV[dt][eid] = slv
        self.slipW[dt][eid] = slp
