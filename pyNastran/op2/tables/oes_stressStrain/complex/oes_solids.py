#pylint: disable=C0301,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import zip, range

from numpy import zeros, concatenate
#from numpy.linalg import eig

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeImagFloats13E


class ComplexSolidArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        #self.cid = {}  # gridGauss

        if is_sort1:
            #sort1
            pass
        else:
            raise NotImplementedError('SORT2')

    def combine(self, results):
        #print('ComplexSolid combine')
        #print('data.shape1 =', self.data.shape)
        #self.data = vstack(data)

        data = [self.nelements] + [result.nelements for result in results]
        self.nelements = sum(data)

        data = [self.ntotal] + [result.ntotal for result in results]
        self.ntotal = sum(data)

        data = [self.element_types3] + [result.element_types3 for result in results]
        self.element_types3 = concatenate(data, axis=0)

        data = [self.element_node] + [result.element_node for result in results]
        self.element_node = concatenate(data, axis=0)

        data = [self.element_cid] + [result.element_cid for result in results]
        self.element_cid = concatenate(data, axis=0)

        data = [self.data] + [result.data for result in results]
        self.data = concatenate(data, axis=1)

        self.data = concatenate(data, axis=1)

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_nnodes(self):
        if self.element_type == 39: # CTETRA
            nnodes = 5
        elif self.element_type == 68: # CPENTA
            nnodes = 7
        elif self.element_type == 67: # CHEXA
            nnodes = 9
        else:
            raise NotImplementedError(self.element_name)
        return nnodes

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = self.get_nnodes()

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.element_types2 = array(self.nelements, dtype='|S8')
        self.element_types3 = zeros((self.nelements, 2), dtype='int32')


        #self.ntotal = self.nelements * nnodes

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((self.ntotal, 2), 'int32')
        self.element_cid = zeros((self.nelements, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)

        if self.result_flag == 0:
            # [oxx, oyy, ozz, txy, tyz, txz]
            self.data = zeros((self.ntimes, self.ntotal, 6), 'complex64')
        else:
            # oxx
            self.data = zeros((self.ntimes, self.ntotal, 1), 'complex64')

    def add_eid_sort1(self, element_num, element_type, dt, eid, cid, ctype, nodef):
        self._times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #self.element_types2[self.ielement] = string_(element_type)   # TODO: save this...
        #self.element_types2[self.ielement] = element_type

        #try:
        if self.ielement < self.nelements:
            self.element_cid[self.ielement] = [eid, cid]
            self.element_types3[self.ielement, :] = [element_num, nodef]
        #except IndexError:
            #pass
            #print('element_types3', self.element_types3)

        #self.node_element_cid[self.itotal] = []
        #self.element_node[self.itotal, :] = [eid, 0]  # 0 is center
        #print("etype=%s ctype=%s nodef=%s" % (element_type, ctype, nodef))
        self.ielement += 1
        #self.itotal += 1

    def add_node_sort1(self, dt, eid, grid, inode, ex, ey, ez, etxy, etyz, etzx):
        #ielem = self.ielement - 1
        #print('eid=%i grid=%i exx=%s' % (eid, grid, str(ex)))
        #print('data[%s, %s, :] = %s dt=%s' % (self.itime, self.itotal, ex, dt))
        #if dt != 0.0:
            #asfd
        #print('data[%s, %s, :] = %s' % (self.itime, self.itotal, [ex, ey, ez, etxy, etyz, etzx]))
        #print('element_node[%s, %s] = [%s, %s]' % (ielem, inode, eid, grid))

        if self.result_flag == 0:
            self.data[self.itime, self.itotal, :] = [ex, ey, ez, etxy, etyz, etzx]
        else:
            self.data[self.itime, self.itotal, 0] = ex
        self.element_node[self.itotal, :] = [eid, grid]
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
        #ntotal = self.ntotal
        nnodes = self.element_node.shape[0]
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)

        # write the f06
        ntimes = self.data.shape[0]

        cid = 0
        for itime in range(ntimes):
            dt = self._times[itime]

            #print('eids=', eids)

            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write('\n'.join(msg))

            # TODO: can I get this without a reshape?
            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            ozz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]

            eids2 = self.element_node[:, 0]
            nodes = self.element_node[:, 1]

            # loop over all the elements and nodes
            for deid, node, doxx, doyy, dozz, dtxy, dtyz, dtxz in zip(eids2, nodes, oxx, oyy, ozz, txy, tyz, txz):
                # TODO: cid not supported
                ([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                  oxxi, oyyi, ozzi, txyi, tyzi, txzi,], is_all_zeros) = writeImagFloats13E([doxx, doyy, dozz,
                                                                                            dtxy, dtyz, dtxz], is_mag_phase)
                if node == 0:  # CENTER
                    f.write('0 %12i %11sGRID CS %2i GP\n' % (deid, cid, nnodes))
                    f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('CENTER', oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                    f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('', oxxi, oyyi, ozzi, txyi, tyzi, txzi))
                else:
                    f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (node, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                    f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('', oxxi, oyyi, ozzi, txyi, tyzi, txzi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexSolidStressArray(ComplexSolidArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz']
        return headers

def _get_msgs(self, is_mag_phase, is_sort1):
    if is_mag_phase:
        mag_phase = '                                                          (MAGNITUDE/PHASE)'
    else:
        mag_phase = '                                                          (REAL/IMAGINARY)'

    if self.is_stress():
        base_msg = [
            mag_phase,
            '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------',
            '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
            '', ]
        tetra_msg = ['                 C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )', ]
        hexa_msg = ['                 C O M P L E X   S T R E S S E S   I N   H E X A H E D R O N   E L E M E N T S   ( C H E X A )', ]
        penta_msg = ['                 C O M P L E X   S T R E S S E S   I N   P E N T A H E D R O N   E L E M E N T S   ( C P E N T A )', ]
    else:
        base_msg = [
            mag_phase,
            '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
            '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
            '', ]
        tetra_msg = ['                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',]
        hexa_msg = ['                 C O M P L E X     S T R A I N S   I N   H E X A H E D R O N   E L E M E N T S   ( C H E X A )',]
        penta_msg = ['                 C O M P L E X     S T R A I N S   I N   P E N T A H E D R O N   E L E M E N T S   ( C P E N T A )',]

    tetra_msg += base_msg
    penta_msg += base_msg
    hexa_msg += base_msg
    return tetra_msg, penta_msg, hexa_msg


def get_f06_header(self, is_mag_phase=True, is_sort1=True):
    tetra_msg, penta_msg, hexa_msg = _get_msgs(self, is_mag_phase, is_sort1)

    if self.element_type == 39:  # CTETRA
        return tetra_msg, 4
    elif self.element_type == 67:  # CHEXA
        return hexa_msg, 8
    elif self.element_type == 68:  # CPENTA
        return penta_msg, 6
    else:
        raise NotImplementedError('complex solid stress/strain name=%r Type=%s' % (self.element_name, self.element_type))


class ComplexSolidStrainArray(ComplexSolidArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        return headers




#class ComplexSolidStress(StressObject):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StressObject.__init__(self, data_code, isubcase)

        #self.result_flag = 0
        #self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.cid = {}  # gridGauss
        #if self.result_flag == 0:
            #self.oxx = {}
            #self.oyy = {}
            #self.ozz = {}
            #self.txy = {}
            #self.tyz = {}
            #self.txz = {}
        #elif self.result_flag == 1:
            #self.wmax = {}
            #self.oxx = {}  # dummy

        #if is_sort1:
            #pass
            ##if dt is not None:
                ##self.add = self.add_sort1
                ##self.add_new_eid = self.add_new_eid_sort1
        #else:
            #assert dt is not None
            #raise NotImplementedError('SORT2')
            ##self.add = self.addSort2
            ##self.add_new_eid = self.add_new_eid_sort2

    #def get_stats(self):
        #nelements = len(self.eType)

        #msg = []
        #if self.nonlinear_factor is not None:  # transient
            #ntimes = len(self.oxx)
            #msg.append('  type=%s ntimes=%i nelements=%i\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        #msg.append('  eType, cid, oxx, oyy, ozz, txy, tyz, txz\n  ')
        #msg.append('%s\n  ' % self.element_name)
        #msg += self.get_data_code()
        #return msg

    #def add_f06_data(self, data, transient):
        #if transient is None:
            #if not hasattr(self, 'data'):
                #self.data = []
            #self.data += data
        #else:
            #dt = transient[1]
            #if not hasattr(self, 'data'):
                #self.data = {}
            ##print(self.data)
            #if dt not in self.data:
                #self.data[dt] = []
            #for line in data:
                #self.data[dt] += data

    #def processF06Data(self):
        #raise NotImplementedError()

    #def delete_transient(self, dt):
        #del self.oxx[dt]
        #del self.oyy[dt]
        #del self.ozz[dt]
        #del self.txy[dt]
        #del self.tyz[dt]
        #del self.txz[dt]

    #def get_transients(self):
        #k = self.oxx.keys()
        #k.sort()
        #return k

    #def add_new_transient(self, dt):
        #"""
        #initializes the transient variables
        #"""
        #if self.result_flag == 0:
            #self.oxx[dt] = {}
            #self.oyy[dt] = {}
            #self.ozz[dt] = {}
            #self.txy[dt] = {}
            #self.tyz[dt] = {}
            #self.txz[dt] = {}
        #elif self.result_flag == 1:
            #self.wmax[dt] = {}
        #else:
            #raise RuntimeError(self.result_flag)


    #def add_eid_sort1(self, element_num, eType, dt, eid, cid, ctype, nodef):
        #assert cid >= -1, cid
        #assert eid >= 0, eid

        ##print "dt=%s eid=%s eType=%s" %(dt,eid,eType)
        #self.eType[eid] = eType
        #self.cid[eid] = cid

        #if self.result_flag == 0:
            #if dt not in self.oxx:
                #self.add_new_transient(dt)
            #assert eid not in self.oxx[dt], self.oxx[dt]
            #self.oxx[dt][eid] = {}
            #self.oyy[dt][eid] = {}
            #self.ozz[dt][eid] = {}
            #self.txy[dt][eid] = {}
            #self.tyz[dt][eid] = {}
            #self.txz[dt][eid] = {}
        #else:
            #if dt not in self.wmax:
                #self.add_new_transient(dt)
            #assert eid not in self.wmax[dt], self.wmax[dt]
            #self.wmax[dt][eid] = {}
        ##msg = "*eid=%s node_id=%s vm=%g" % (eid, node_id, ovm)

    #def add_node_sort1(self, dt, eid, node_id, inode, oxx, oyy, ozz, txy, tyz, tzx):
        ##msg = "eid=%s node_id=%s vm=%g" % (eid, node_id, ovm)
        ##print msg

        #assert isinstance(node_id, int), node_id

        #if self.result_flag == 1:
            #self.wmax[dt][eid][node_id] = oxx
            #return

            #tensor = zeros((3, 3), dtype='complex64')
            #tensor[0, 0] = oxx
            #tensor[1, 1] = oyy
            #tensor[2, 2] = ozz
            #tensor[0, 1] = txy
            #tensor[1, 2] = tyz
            #tensor[0, 2] = tzx
            #tensor[1, 0] = tensor[0, 1]
            #tensor[2, 1] = tensor[1, 2]
            #tensor[2, 0] = tensor[0, 2]

            #wr, v = eig(tensor.real)
            #wi, v = eig(tensor.imag)
            #self.wmax[dt][eid][node_id] = abs(wr + 1j * wi).max()   #close...
            #return

            ##w, v = eig(tensor)
            ##self.wmax[dt][eid][node_id] = abs(w).max()
            ##return

            ##w, v = eig(abs(tensor))  #  close for freq=343
            ##self.wmax[dt][eid][node_id] = abs(w).max()
            ##return

            ##w, v = eig(abs(tensor.real))
            ##self.wmax[dt][eid][node_id] = w.max()
            ##return
            ##w, v = eig(tensor)
            ##abs_w = abs(w)
            ##iw = list(abs_w).index(abs_w.max())
            ##self.wmax[dt][eid][node_id] = w[iw]
            ##return

            ##wreal = eigvalsh(tensor.real)
            ##wimag = eigvalsh(tensor.imag)
            ##wmaxr = 0.
            ##wmaxi = 0.
            ##for ii in range(3):
                ##if abs(wreal[ii]) > abs(wmaxr):
                    ##wmaxr = wreal[ii]
                ##if abs(wimag[ii]) > abs(wmaxi):
                    ##wmaxi = wimag[ii]
            ##wmax = wmaxr + 1j * wmaxi
            ##self.wmax[dt][eid][node_id] = wmax
        #else:
            ##print("eid=%s nid=%s oyy=%s" %(eid, node_id, oyy))
            #self.oxx[dt][eid][node_id] = oxx
            #self.oyy[dt][eid][node_id] = oyy
            #self.ozz[dt][eid][node_id] = ozz

            #self.txy[dt][eid][node_id] = txy
            #self.tyz[dt][eid][node_id] = tyz
            #self.txz[dt][eid][node_id] = tzx

    #def get_headers(self):
        #headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz']
        #return headers

    #def directional_vectors(self, oxx, oyy, ozz, txy, tyz, txz):
        #A = [[oxx, txy, txz],
             #[txy, oyy, tyz],
             #[txz, tyz, ozz]]
        #(Lambda, v) = eig(A)  # we can't use a hermitian matrix
        #return v

    #def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #tetra_msg, penta_msg, hexa_msg = _get_msgs(self, is_mag_phase, is_sort1)
        #dts = list(self.oxx.keys())
        #dts.sort()
        #dt0 = dts[0]
        #eids = sorted(self.oxx[dt0].keys())

        #if self.element_type == 39:  # CTETRA
            #for dt in dts:
                #self._write_element_transient('CTETRA', 4, eids, dt, header, tetra_msg, f, is_mag_phase)
                #f.write(page_stamp % page_num)
                #page_num += 1
        #elif self.element_type == 67:  # CHEXA
            #for dt in dts:
                #self._write_element_transient('CHEXA', 8, eids, dt, header, hexa_msg, f, is_mag_phase)
                #f.write(page_stamp % page_num)
                #page_num += 1
        #elif self.element_type == 68:  # CPENTA
            #for dt in dts:
                #self._write_element_transient('CPENTA', 6, eids, dt, header, penta_msg, f, is_mag_phase)
                #f.write(page_stamp % page_num)
                #page_num += 1
        #else:
            #raise NotImplementedError('complex solid stress/strain name=%r Type=%s' % (self.element_name, self.element_type))
        #return page_num - 1

    #def _write_element_transient(self, element_name, nnodes, eids, dt, header, msg, f, is_mag_phase):
        #dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
        #header[1] = dt_line
        #msg = header + msg

        #cid = 0
        #f.write('\n'.join(msg))
        #for eid in eids:
            #node_ids = sorted(self.oxx[dt][eid].keys())
            #f.write('0 %12i %11sGRID CS %2i GP\n' % (eid, cid, nnodes))
            #for inode in node_ids:
                #oxx = self.oxx[dt][eid][inode]
                #oyy = self.oyy[dt][eid][inode]
                #ozz = self.ozz[dt][eid][inode]
                #txy = self.txy[dt][eid][inode]
                #tyz = self.tyz[dt][eid][inode]
                #txz = self.txz[dt][eid][inode]
                #([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                  #oxxi, oyyi, ozzi, txyi, tyzi, txzi,], is_all_zeros) = writeImagFloats13E([oxx, oyy, ozz,
                                                                                            #txy, tyz, txz], is_mag_phase)

                #f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (inode, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                #f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('', oxxi, oyyi, ozzi, txyi, tyzi, txzi))


#class ComplexSolidStrain(StrainObject):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #StrainObject.__init__(self, data_code, isubcase)

        #self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.cid = {}  # gridGauss
        #self.exx = {}
        #self.eyy = {}
        #self.ezz = {}
        #self.exy = {}
        #self.eyz = {}
        #self.exz = {}

        ##self.dt = dt
        #if is_sort1:
            #if dt is not None:
                #pass
        #else:
            #assert dt is not None
            #raise NotImplementedError('SORT2')

    #def get_stats(self):
        #nelements = len(self.eType)

        #msg = []
        #if self.nonlinear_factor is not None:  # transient
            #ntimes = len(self.exx)
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__, nelements))
        #msg.append('  eType, cid, exx, eyy, ezz, exy, eyz, exz\n  ')
        #msg.append('%s\n  ' % self.element_name)
        #msg = msg + self.get_data_code()
        #return msg

    #def add_f06_data(self, data, transient):
        #if transient is None:
            #if not hasattr(self, 'data'):
                #self.data = []
            #self.data += data
        #else:
            #dt = transient[1]
            #if not hasattr(self, 'data'):
                #self.data = {}
            ##print(self.data)
            #if dt not in self.data:
                #self.data[dt] = []
            #for line in data:
                #self.data[dt] += data

    #def processF06Data(self):
        #raise NotImplementedError()

    #def delete_transient(self, dt):
        #del self.exx[dt]
        #del self.eyy[dt]
        #del self.ezz[dt]
        #del self.exy[dt]
        #del self.eyz[dt]
        #del self.exz[dt]

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
        #self.ezz[dt] = {}
        #self.exy[dt] = {}
        #self.eyz[dt] = {}
        #self.exz[dt] = {}

    #def add_eid_sort1(self, element_num, eType, dt, eid, cid, ctype, nodef):
        #assert cid >= -1, cid
        #assert eid >= 0, eid

        #if dt not in self.exx:
            #self.add_new_transient(dt)
        #assert eid not in self.exx[dt], self.exx[dt]

        #self.eType[eid] = eType
        #self.cid[eid] = cid
        #self.exx[dt][eid] = {}
        #self.eyy[dt][eid] = {}
        #self.ezz[dt][eid] = {}
        #self.exy[dt][eid] = {}
        #self.eyz[dt][eid] = {}
        #self.exz[dt][eid] = {}

        ##msg = "*eid=%s node_id=%s vm=%g" % (eid, node_id, ovm)
        ##print msg

    #def add_node_sort1(self, dt, eid, node_id, inode, ex, ey, ez, etxy, etyz, etzx):
        ##msg = "eid=%s node_id=%s vm=%g" % (eid, node_id, evm)
        #assert isinstance(node_id, int), node_id
        ##print "eid=%s nid=%s exx=%s" %(eid, node_id, exx)
        #self.exx[dt][eid][node_id] = ex
        #self.eyy[dt][eid][node_id] = ey
        #self.ezz[dt][eid][node_id] = ez

        #self.exy[dt][eid][node_id] = etxy
        #self.eyz[dt][eid][node_id] = etyz
        #self.exz[dt][eid][node_id] = etzx

    #def get_headers(self):
        #headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        #return headers

    #def directional_vectors(self, exx, eyy, ezz, exy, eyz, exz):
        #A = [[exx, exy, exz],
             #[exy, eyy, eyz],
             #[exz, eyz, ezz]]
        #(Lambda, v) = eig(A)  # we can't use a hermitian matrix
        #return v

    #def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #tetra_msg, hexa_msg, penta_msg = _get_msgs(self, is_mag_phase, is_sort1)

        #dts = list(self.exx.keys())
        #dts.sort()
        #dt0 = dts[0]
        #eids = sorted(self.exx[dt0].keys())

        #for dt in dts:
            #if self.element_type == 39:  # CTETRA
                #self._write_element_transient('CTETRA', 4, eids, dt, header, tetra_msg, f, is_mag_phase)
                #f.write(page_stamp % page_num)
                #page_num += 1
            #elif self.element_type == 68:  # CPENTA
                #self._write_element_transient('CPENTA', 6, eids, dt, header, penta_msg, f, is_mag_phase)
                #f.write(page_stamp % page_num)
                #page_num += 1
            #elif self.element_type == 67:  # CHEXA
                #self._write_element_transient('CHEXA', 8, eids, dt, header, hexa_msg, f, is_mag_phase)
                #f.write(page_stamp % page_num)
                #page_num += 1
            #else:
                #raise NotImplementedError('complex solid stress/strain name=%r Type=%s' % (self.element_name, self.element_type))
        #return page_num - 1

    #def _write_element_transient(self, element_name, nnodes, eids, dt, header, msg, f, is_mag_phase):
        #dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
        #header[1] = dt_line
        #msg = header + msg

        #f.write('\n'.join(msg))
        #cid = 0
        #for eid in eids:
            #node_ids = sorted(self.exx[dt][eid].keys())

            #f.write('0 %12i %11sGRID CS %2i GP\n' % (eid, cid, nnodes))  ## TODO: cid
            #for inode in node_ids:
                #oxx = self.exx[dt][eid][inode]
                #oyy = self.eyy[dt][eid][inode]
                #ozz = self.ezz[dt][eid][inode]
                #txy = self.exy[dt][eid][inode]
                #tyz = self.eyz[dt][eid][inode]
                #txz = self.exz[dt][eid][inode]
                #([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                  #oxxi, oyyi, ozzi, txyi, tyzi, txzi,], is_all_zeros) = writeImagFloats13E([oxx, oyy, ozz,
                                                                                            #txy, tyz, txz], is_mag_phase)

                #f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (inode, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                #f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('', oxxi, oyyi, ozzi, txyi, tyzi, txzi))
