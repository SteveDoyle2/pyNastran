from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from six.moves import zip, range
from numpy import zeros, searchsorted

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import _eigenvalue_header, get_key0

class RealShearArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]
        self.nelements = 0  # result specific

        if is_sort1:
            self.add_new_eid = self.add_new_eid_sort1
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
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self.times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        # [max_shear, avg_shear, margin]
        self.data = zeros((self.ntimes, self.ntotal, 3), dtype='float32')

    def add_new_eid_sort1(self, dt, eid, max_shear, avg_shear, margin):
        """
        ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY
          ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN
            328        1.721350E+03   1.570314E+03   7.2E+01
        """
        self.times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [max_shear, avg_shear, margin]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self):
        raise NotImplementedError('CSHEAR...')

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element_node[:, 0])
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg_temp = self.get_f06_header()

        # write the f06
        (ntimes, ntotal, three) = self.data.shape

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        for itime in range(ntimes):
            dt = self.times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            max_shear = self.data[itime, :, 0]
            avg_shear = self.data[itime, :, 1]
            margin = self.data[itime, :, 2]

            # loop over all the elements
            out = []
            for eid, max_sheari, avg_sheari, margini in zip(eids, max_shear, avg_shear, margin):
                #([max_sheari, avg_sheari, margini],
                #is_all_zeros) = writeFloats13E([max_sheari, avg_sheari, margini])
                out.append([eid, max_sheari, avg_sheari, margini])

            for i in range(0, nwrite, 2):
                out_line = '      %8i   %13s  %10.4E %13s  %8i   %13s  %10.4E %s\n' % (tuple(out[i] + out[i + 1]))
                f.write(out_line)
            if is_odd:
                out_line = '      %8i   %13s  %10.4E %s\n' % tuple(out[-1])
                f.write(out_line)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealShearStressArray(RealShearArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['max_shear', 'avg_shear', 'margin']
        return headers

    def get_f06_header(self):
        msg = [
            '                                     S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )\n'
            '      ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY\n'
            '        ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN\n'
           #'          328        1.721350E+03   1.570314E+03   7.2E+01'
            ]
        return msg


class RealShearStrainArray(RealShearArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['max_shear', 'avg_shear', 'margin']
        return headers

    def get_f06_header(self):
        msg = [
            '                                     S T R A I N S   I N   S H E A R   P A N E L S      ( C S H E A R )\n'
            '      ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY\n'
            '        ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN\n'
           #'          328        1.721350E+03   1.570314E+03   7.2E+01'
            ]
        return msg


class RealShearStress(StressObject):
    """
    ::

      # format_code=1 sort_code=0 stressCode=0
                                     S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )
      ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY
        ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN
          328        1.721350E+03   1.570314E+03   7.2E+01
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CSHEAR'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}

        self.getLength = self.getLength
        self.isImaginary = False
        #if dt is not None:
        #    self.add_new_transient = self.add_new_transient
        #    self.add_new_eid       = self.addNewEidTransient
        #else:
        #    self.add_new_eid = self.add_new_eid

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.maxShear)
            s0 = get_key0(self.maxShear)
            nelements = len(self.maxShear[s0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.maxShear)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, maxShear, avgShear, margin\n')
        return msg

    def delete_transient(self, dt):
        del self.maxShear[dt]
        del self.avgShear[dt]
        del self.margin[dt]

    def get_transients(self):
        k = self.maxShear.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.maxShear[dt] = {}
        self.avgShear[dt] = {}
        self.margin[dt] = {}

    def getLength(self):
        return (16, 'fff')

    def add_f06_data(self, data, dt):
        if dt:
            if dt not in self.maxShear:
                self.maxShear[dt] = {}
                self.avgShear[dt] = {}
                self.margin[dt] = {}
            for datai in data:
                (eid, max_shear, avg_shear, margin) = datai
                self.maxShear[dt][eid] = max_shear
                self.avgShear[dt][eid] = avg_shear
                self.margin[dt][eid] = margin
                return

        for datai in data:
            (eid, max_shear, avg_shear, margin) = datai
            self.maxShear[eid] = max_shear
            self.avgShear[eid] = avg_shear
            self.margin[eid] = margin

    def add_new_eid(self, dt, eid, maxShear, avgShear, margin):
        #print "Rod Stress add..."
        assert isinstance(eid, int)
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}
        self.maxShear[eid] = maxShear
        self.avgShear[eid] = avgShear
        self.margin[eid] = margin

    def add_new_eid_sort1(self, dt, eid, maxShear, avgShear, margin):
        if dt not in self.maxShear:
            self.add_new_transient(dt)
        assert isinstance(eid, int)
        assert eid >= 0, eid
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin


class RealShearStrain(StrainObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = 'CSHEAR'
        #raise Exception('not supported...CSHEAR strain')
        self.code = [self.format_code, self.sort_code, self.s_code]
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.maxShear)
            s0 = get_key0(self.maxShear)
            nelements = len(self.maxShear[s0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.maxShear)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, maxShear, avgShear, margin\n')
        return msg

    def delete_transient(self, dt):
        del self.maxShear[dt]
        del self.avgShear[dt]
        del self.margin[dt]

    def get_transients(self):
        k = self.maxShear.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        .. note:: make sure you set self.dt first
        """
        self.dt = dt
        self.maxShear[dt] = {}
        self.avgShear[dt] = {}
        self.margin[dt] = {}

    def add_f06_data(self, data, dt):
        if dt:
            if dt not in self.maxShear:
                self.maxShear[dt] = {}
                self.avgShear[dt] = {}
                self.margin[dt] = {}
            for datai in data:
                (eid, max_shear, avg_shear, margin) = datai
                self.maxShear[dt][eid] = max_shear
                self.avgShear[dt][eid] = avg_shear
                self.margin[dt][eid] = margin
                return

        for datai in data:
            (eid, max_shear, avg_shear, margin) = datai
            self.maxShear[eid] = max_shear
            self.avgShear[eid] = avg_shear
            self.margin[eid] = margin

    #def add_new_eid(self, dt, eid, axial, SMa, torsion, SMt):
        #raise NotImplementedError()
        ##(axial, SMa, torsion, SMt) = out
        #print "Rod Strain add..."
        #assert eid >= 0, eid
        #self.eType = self.eType
        #self.maxShear[eid] = axial
        #self.avgShear[eid] = SMa
        #self.margin[eid] = torsion

    def add_new_eid_sort1(self, dt, eid, maxShear, avgShear, margin):
        #(maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.add_new_transient(dt)
        assert eid >= 0, eid

        #self.eType[eid] = self.element_type
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin

    def add_new_eid_sort2(self, eid, dt, maxShear, avgShear, margin):
        #(maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.add_new_transient(dt)
        assert eid >= 0, eid

        #self.eType[eid] = self.element_type
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin