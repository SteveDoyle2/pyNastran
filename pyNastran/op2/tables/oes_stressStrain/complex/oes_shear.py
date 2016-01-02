from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import get_key0


class ComplexShearArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_nnodes(self):
        return get_nnodes(self)

    def build(self):
        #print('data_code = %s' % self.data_code)
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.ntotal = self.nelements * nnodes

        self.element = zeros((self.nelements, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)

        # [max_shear, avg_shear]
        self.data = zeros((self.ntimes, self.ntotal, 2), 'complex64')

    def add_sort1(self, dt, eid, max_shear, avg_shear):
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal] = [max_shear, avg_shear]
        self.element[self.itotal, :] = eid
        #self.ielement += 1
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
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  data: [ntimes, nnodes, 2] where 2=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        raise NotImplementedError('ComplexShearArray')
        msg_temp, nnodes, is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1)

        ntimes = self.data.shape[0]
        eids = self.element
        if self.is_sort1():
            if is_sort1:
                for itime in range(ntimes):
                    dt = self._times[itime]

                    dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                    header[1] = dtLine
                    msg = header + msg_temp
                    f.write('\n'.join(msg))

                    max_shear = self.data[itime, :, 0]
                    avg_shear = self.data[itime, :, 1]
                    for eid, max_sheari, avg_sheari in zip(eids, max_shear, avg_shear):
                        ([rmax_shear, imax_shear, ravg_shear, iavg_shear
                          ,], is_all_zeros) = writeImagFloats13E([max_sheari, avg_sheari], is_magnitude_phase)

                        f.write('   %6s   %-13s / %-13s     %-13s / %-13s\n' % (
                            eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                    f.write(page_stamp % page_num)
                    page_num += 1
            else:
                times = self._times
                for ieid, eid in enumerate(eids):
                    max_shear = self.data[:, ieid, 0].ravel()
                    avg_shear = self.data[:, ieid, 1].ravel()
                    for itime, max_sheari, avg_sheari in zip(times, max_shear, avg_shear):
                        ([rmax_shear, imax_shear, ravg_shear, iavg_shear
                          ,], is_all_zeros) = writeImagFloats13E([max_sheari, avg_sheari], is_magnitude_phase)

                        f.write('   %6s   %-13s / %-13s     %-13s / %-13s\n' % (
                            eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                    f.write(page_stamp % page_num)
                    page_num += 1
        else:
            raise NotImplementedError('ComplexShearArray-sort2')
        return page_num - 1


class ComplexShearStressArray(ComplexShearArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        return ['max_shear', 'avg_shear']

class ComplexShearStrainArray(ComplexShearArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain(), self.stress_bits

    def _get_headers(self):
        return ['max_shear', 'avg_shear']

#class ComplexShearStress(StressObject):
    #"""
    #::

      ## format_code=1 sort_code=0 stressCode=0
                                     #S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )
      #ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY
        #ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN
          #328        1.721350E+03   1.570314E+03   7.2E+01
    #"""
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StressObject.__init__(self, data_code, isubcase)
        #self.eType = 'CSHEAR'

        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.maxShear = {}
        #self.avgShear = {}

        #self.isImaginary = False

        #self.dt = dt
        #if is_sort1:
            #if dt is not None:
                ##self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        #else:
            #assert dt is not None
            ##self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2

    #def get_stats(self):
        #msg = self.get_data_code()
        #if self.dt is not None:  # transient
            #ntimes = len(self.maxShear)
            #s0 = get_key0(self.maxShear)
            #nelements = len(self.maxShear[s0])
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #nelements = len(self.maxShear)
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
        #msg.append('  eType, maxShear, avgShear\n')
        #return msg

    #def delete_transient(self, dt):
        #del self.maxShear[dt]
        #del self.avgShear[dt]

    #def get_transients(self):
        #k = self.maxShear.keys()
        #k.sort()
        #return k

    #def add_new_transient(self, dt):
        #"""
        #initializes the transient variables
        #"""
        #self.dt = dt
        #self.maxShear[dt] = {}
        #self.avgShear[dt] = {}

    #def add_f06_data(self, data, dt):
        #if dt:
            #if dt not in self.maxShear:
                #self.maxShear[dt] = {}
                #self.avgShear[dt] = {}
            #for datai in data:
                #(eid, max_shear, avg_shear) = datai
                #self.maxShear[dt][eid] = max_shear
                #self.avgShear[dt][eid] = avg_shear
                #return

        #for datai in data:
            #(eid, max_shear, avg_shear) = datai
            #self.maxShear[eid] = max_shear
            #self.avgShear[eid] = avg_shear

    #def add_new_eid_sort1(self, dt, eid, max_shear, avg_shear):
        #if dt not in self.maxShear:
            #self.add_new_transient(dt)
        #assert isinstance(eid, int)
        #assert eid >= 0, eid
        #self.maxShear[dt][eid] = max_shear
        #self.avgShear[dt][eid] = avg_shear


#class ComplexShearStrain(StrainObject):

    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StrainObject.__init__(self, data_code, isubcase)
        #self.eType = 'CSHEAR'
        ##raise Exception('not supported...CSHEAR strain')
        #self.code = [self.format_code, self.sort_code, self.s_code]
        #self.maxShear = {}
        #self.avgShear = {}

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
        #if self.dt is not None:  # transient
            #ntimes = len(self.maxShear)
            #s0 = get_key0(self.maxShear)
            #nelements = len(self.maxShear[s0])
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #nelements = len(self.maxShear)
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
        #msg.append('  eType, maxShear, avgShear\n')
        #return msg

    #def delete_transient(self, dt):
        #del self.maxShear[dt]
        #del self.avgShear[dt]

    #def get_transients(self):
        #k = self.maxShear.keys()
        #k.sort()
        #return k

    #def add_new_transient(self, dt):
        #"""
        #initializes the transient variables
        #.. note:: make sure you set self.dt first
        #"""
        #self.dt = dt
        #self.maxShear[dt] = {}
        #self.avgShear[dt] = {}

    #def add_f06_data(self, data, dt):
        #if dt:
            #if dt not in self.maxShear:
                #self.maxShear[dt] = {}
                #self.avgShear[dt] = {}
            #for datai in data:
                #(eid, max_shear, avg_shear) = datai
                #self.maxShear[dt][eid] = max_shear
                #self.avgShear[dt][eid] = avg_shear
                #return

        #for datai in data:
            #(eid, max_shear, avg_shear) = datai
            #self.maxShear[eid] = max_shear
            #self.avgShear[eid] = avg_shear

    #def add_new_eid_sort1(self, dt, eid, max_shear, avg_shear):
        #if dt not in self.maxShear:
            #self.add_new_transient(dt)
        #assert eid >= 0, eid

        ##self.eType[eid] = self.element_type
        #self.maxShear[dt][eid] = max_shear
        #self.avgShear[dt][eid] = avg_shear

    #def add_new_eid_sort2(self, eid, dt, max_shear, avg_shear):
        #if dt not in self.maxShear:
            #self.add_new_transient(dt)
        #assert eid >= 0, eid

        ##self.eType[eid] = self.element_type
        #self.maxShear[dt][eid] = max_shear
        #self.avgShear[dt][eid] = avg_shear
