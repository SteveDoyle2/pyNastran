#pylint: disable=C0301,C0111,R0921
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip
from collections import defaultdict

from numpy import zeros, array, string_, searchsorted, where, argwhere, ravel
from numpy.linalg import eig

from ..real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeImagFloats13E


class ComplexSolid(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        #self.cid = {}  # gridGauss

        if is_sort1:
            #sort1
            pass
        else:
            raise NotImplementedError('SORT2')

    def _get_msgs(self, is_mag_phase):
        raise NotImplementedError()

    def build(self):
        if self.is_built:
            self.itotal = 0
            self.ielement = 0
            return

        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self.times = zeros(self.ntimes, 'float32')
        #self.element_types2 = array(self.nelements, dtype='|S8')
        self.element_types3 = zeros((self.nelements, 2), dtype='int32')

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((self.ntotal, 2), 'int32')
        self.element_cid = zeros((self.nelements, 2), 'int32')

        #[oxx, oyy, ozz, txy, tyz, txz]
        self.data = zeros((self.ntimes, self.ntotal, 6), 'complex64')

    def add_eid_sort1(self, element_num, element_type, dt, eid, cid, ctype, nodef):
        self.times[self.itime] = dt
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
        self.element_node[self.itotal, :] = [eid, 0]  # 0 is center
        #print("etype=%s ctype=%s nodef=%s" % (element_type, ctype, nodef))
        self.ielement += 1
        #self.itotal += 1
        #self.data

    def add_node_sort1(self, dt, eid, grid, ex, ey, ez, etxy, etyz, etzx):
        #print('eid=%i grid=%i exx=%s' % (eid, grid, str(ex)))
        self.data[self.itime, self.itotal, :] = [ex, ey, ez, etxy, etyz, etzx]
        #print('data[%s, %s, :] = %s' % (self.itime, self.itotal, str(self.data[self.itime, self.itotal, :])))
        self.element_node[self.itotal] = [eid, grid]
        self.itotal += 1
        #self.data

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
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
        #msg.append('  data.shape=%s' % str(self.data.shape))
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n  ' % str(', '.join(self.getHeaders())))

        msg.append(', '.join(set(self.eType.values())))
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        tetra_msg, penta_msg, hexa_msg = self._get_msgs(is_mag_phase)

        # TODO: this should use numpy/scipy
        d1 = defaultdict(int)
        d2 = defaultdict(list)

        i = 0
        for (etype, nnodes) in self.element_types3:
            d1[(etype, nnodes)] += 1  # count???
            d2[(etype, nnodes)].append(self.element_cid[i, 0])  # indexs
            i += 1
        #print("d1=%s" % d1)
        #print("d2=%s" % d2)
        elements = {}
        TETRA = {}
        HEXA = {}
        ENTA = {}
        for key in d1:
            (eType, nnodes) = key
            if eType == 39:  # TODO: doesn't support CPENTA, CHEXA
                TETRA['CTETRA'+str(nnodes)] = array(d2[key])
            #elif eType == 39:
                #elements['CTETRA'+str(nnodes)] = d2[key]
            else:
                raise NotImplementedError(eType)
        #print(elements)
        return (tetra_msg, hexa_msg, penta_msg,
                TETRA, HEXA, PENTA)

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([argwhere(self.element_node[:, 0]==eid) for eid in eids])
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, pageNum=1, f=None, is_mag_phase=False):
        (tetra_msg, hexa_msg, penta_msg,
         TETRA, HEXA, PENTA) = self.get_f06_header(is_mag_phase)

        # get the ieids by element type so we don't need to regenerate it inside the dt loop
        ieids_map = {}
        for dictionary in [TETRA, HEXA, PENTA]:
            # don't worry about optimizing this length ~2 dictionary
            for Type, eids in sorted(dictionary.iteritems()):  # CTETRA4, CTETRA10, CTETRA8, etc.
                ieids = self.eid_to_element_node_index(eids)
                ieids_map[Type] = ieids

        # write the f06
        (ntimes, ntotal, six) = self.data.shape
        for itime in xrange(ntimes):
            dt = self.times[itime]  ## TODO: rename this...
            for dictionary, msg_temp in zip([TETRA, HEXA, PENTA], [tetra_msg, hexa_msg, penta_msg,]):
                for Type, eids in sorted(dictionary.iteritems()):
                    #print('eids=', eids)

                    dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                    header[1] = dtLine
                    msg = header + msg_temp
                    f.write('\n'.join(msg))

                    if 'CTETRA' in Type or 'CPENTA' in Type:
                        nnodes = int(Type[6:])  # CTETRA4 / CTETRA10 / CPENTA6 / CPENTA15
                    elif 'CHEXA' in Type:
                        nnodes = int(Type[5:])  # CHEXA8 / CHEXA20
                    else:
                        raise NotImplementedError('complex solid stress/strain Type=%r' % Type)

                    # we know all the CTETRA eids, so let's
                    # get all the CTETRA indicies that we're going to write
                    ieids = ieids_map[Type]  # faster way...
                    #ieids = self.eid_to_element_node_index(eids)

                    #ieids = self.get_element_index(eids)
                    #eid = self.element_node[:, 0]
                    #grid = self.element_node[:, 1]
                    #print('ieids=', ieids)
                    #print('eid=', eid)
                    #print('grid=', grid)
                    #print('**', self.element_node[ieids, 0])
                    #oxx, oyy, ozz, txy, tyz, txz = self.data[itime, ieids, :]
                    n = len(ieids)
                    #print('itime=', itime)

                    # TODO: can I get this without a reshape?
                    oxx = self.data[itime, ieids, 0].reshape(n)
                    oyy = self.data[itime, ieids, 1].reshape(n)
                    ozz = self.data[itime, ieids, 2].reshape(n)
                    txy = self.data[itime, ieids, 3].reshape(n)
                    tyz = self.data[itime, ieids, 4].reshape(n)
                    txz = self.data[itime, ieids, 5].reshape(n)

                    eids2 = self.element_node[ieids, 0]
                    nodes = self.element_node[ieids, 1]
                    #print('eids2 =', eids2)
                    #print('nodes =', nodes)
                    #print('oxx =', oxx)
                    # loop over all the elements and nodes
                    for deid, node, doxx, doyy, dozz, dtxy, dtyz, dtxz in izip(eids2, nodes, oxx, oyy, ozz, txy, tyz, txz):
                        # TODO: cid not supported
                        ([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                          oxxi, oyyi, ozzi, txyi, tyzi, txzi,], isAllZeros) = writeImagFloats13E([doxx, doyy, dozz,
                                                                                                  dtxy, dtyz, dtxz], is_mag_phase)
                        if node == 0:  # CENTER
                            f.write('0 %12i %11sGRID CS %2i GP\n' % (deid, 0, nnodes))
                            f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('CENTER', oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                            f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('',       oxxi, oyyi, ozzi, txyi, tyzi, txzi))
                        else:
                            f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (node, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                            f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('',   oxxi, oyyi, ozzi, txyi, tyzi, txzi))
                    #self.element_types3[ielem]
                if len(dictionary.keys()):
                    f.write(page_stamp % pageNum)
                    pageNum += 1
        return pageNum - 1


class ComplexSolidStressVector(ComplexSolid, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSolid.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def getHeaders(self):
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz']
        return headers

    def _get_msgs(self, is_mag_phase):
        if is_mag_phase:
            base_msg = [
                '                                                          (MAGNITUDE/PHASE)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '', ]

        else:
            base_msg = [
            '                                                          (REAL/IMAGINARY)',
            '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------',
            '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
            '', ]
        tetra_msg = ['                 C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )', ]
        hexa_msg  = ['                 C O M P L E X   S T R E S S E S   I N   H E X A H E D R O N   E L E M E N T S   ( C H E X A )', ]
        penta_msg = ['                 C O M P L E X   S T R E S S E S   I N   P E N T A H E D R O N   E L E M E N T S   ( C P E N T A )', ]
        tetra_msg += base_msg
        penta_msg += base_msg
        hexa_msg += base_msg
        return tetra_msg, penta_msg, hexa_msg


class ComplexSolidStrainVector(ComplexSolid, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSolid.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def getHeaders(self):
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        return headers

    def _get_msgs(self, is_mag_phase):
        if is_mag_phase:
            base_msg = [
                '                                                          (MAGNITUDE/PHASE)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        else:
            base_msg = [
                '                                                          (REAL/IMAGINARY)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        tetra_msg = ['                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',]
        hexa_msg  = ['                 C O M P L E X     S T R A I N S   I N   H E X A H E D R O N   E L E M E N T S   ( C H E X A )',]
        penta_msg = ['                 C O M P L E X     S T R A I N S   I N   P E N T A H E D R O N   E L E M E N T S   ( C P E N T A )',]
        tetra_msg += base_msg
        penta_msg += base_msg
        hexa_msg += base_msg
        return tetra_msg, penta_msg, hexa_msg


class ComplexSolidStressObject(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)

        self.eType = {}
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.cid = {}  # gridGauss
        self.oxx = {}
        self.oyy = {}
        self.ozz = {}
        self.txy = {}
        self.tyz = {}
        self.txz = {}

        if is_sort1:
            pass
            #if dt is not None:
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            raise NotImplementedError('SORT2')
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = []
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid, oxx, oyy, ozz, txy, tyz, txz\n  ')
        msg.append(', '.join(set(self.eType.values()))+'\n  ')
        msg += self.get_data_code()
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            if not hasattr(self, 'data'):
                self.data = []
            self.data += data
        else:
            dt = transient[1]
            if not hasattr(self, 'data'):
                self.data = {}
            #print(self.data)
            if dt not in self.data:
                self.data[dt] = []
            for line in data:
                self.data[dt] += data

    def processF06Data(self):
        raise NotImplementedError()

    def delete_transient(self, dt):
        del self.oxx[dt]
        del self.oyy[dt]
        del self.ozz[dt]
        del self.txy[dt]
        del self.tyz[dt]
        del self.txz[dt]

    def get_transients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.ozz[dt] = {}
        self.txy[dt] = {}
        self.tyz[dt] = {}
        self.txz[dt] = {}

    def add_eid_sort1(self, element_num, eType, dt, eid, cid, ctype, nodef):
        assert cid >= 0
        assert eid >= 0

        #print "dt=%s eid=%s eType=%s" %(dt,eid,eType)
        if dt not in self.oxx:
            self.add_new_transient(dt)
        assert eid not in self.oxx[dt], self.oxx[dt]

        self.eType[eid] = eType
        self.cid[eid] = cid
        self.oxx[dt][eid] = {}
        self.oyy[dt][eid] = {}
        self.ozz[dt][eid] = {}
        self.txy[dt][eid] = {}
        self.tyz[dt][eid] = {}
        self.txz[dt][eid] = {}

        #msg = "*eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
        #print msg
        #if nodeID == 0:
            #raise ValueError(msg)

    def add_node_sort1(self, dt, eid, nodeID, oxx, oyy, ozz, txy, tyz, tzx):
        #msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
        #print msg

        if nodeID == 0:
            nodeID = 'CENTER'
        #print("eid=%s nid=%s oyy=%s" %(eid,nodeID, oyy))
        self.oxx[dt][eid][nodeID] = oxx
        self.oyy[dt][eid][nodeID] = oyy
        self.ozz[dt][eid][nodeID] = ozz

        self.txy[dt][eid][nodeID] = txy
        self.tyz[dt][eid][nodeID] = tyz
        self.txz[dt][eid][nodeID] = tzx

    def getHeaders(self):
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz']
        return headers

    def directionalVectors(self, oxx, oyy, ozz, txy, tyz, txz):
        A = [[oxx, txy, txz],
             [txy, oyy, tyz],
             [txz, tyz, ozz]]
        (Lambda, v) = eig(A)  # we can't use a hermitian matrix
        return v

    def getF06_Header(self, is_mag_phase):
        if is_mag_phase:
            tetra_msg = [
                '                 C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (MAGNITUDE/PHASE)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        else:
            tetra_msg = [
                '                 C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (REAL/IMAGINARY)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        penta_msg = tetra_msg  # TODO: this isnt done
        hexa_msg = tetra_msg  # TODO: this isnt done

        tetra_eids = []
        hexa_eids = []
        penta_eids = []

        tetra10_eids = []
        hexa20_eids = []
        penta15_eids = []
        for eid, eType in sorted(self.eType.iteritems()):
            if eType in ['CTETRA4']:
                tetra_eids.append(eid)
            elif eType in ['CTETRA10']:
                tetra10_eids.append(eid)

            elif eType in ['CPENTA']:
                penta_eids.append(eid)
            elif eType in ['CPENTA15']:
                penta_eids.append(eid)

            elif eType in ['CHEXA8']:
                hexa_eids.append(eid)
            elif eType in ['CHEXA20']:
                hexa20_eids.append(eid)

            else:
                raise NotImplementedError('eType=%r' % eType)
        return (tetra_msg, hexa_msg, penta_msg,
                tetra_eids, hexa_eids, penta_eids,
                tetra10_eids, hexa20_eids, penta15_eids)

    def write_f06(self, header, page_stamp, pageNum=1, f=None, is_mag_phase=False):
        (tetra_msg, hexa_msg, penta_msg,
         tetra_eids, hexa_eids, penta_eids,
         tetra10_eids, hexa20_eids, penta15_eids) = self.getF06_Header(is_mag_phase)
        dts = self.oxx.keys()
        for dt in sorted(dts):
            if tetra_eids:
                self.write_element_transient('CTETRA', 4, tetra_eids, dt, header, tetra_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if tetra10_eids:
                self.write_element_transient('CTETRA', 10, tetra10_eids, dt, header, tetra_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

            if hexa_eids:
                self.write_element_transient('CHEXA',  8,  hexa_eids, dt, header, hexa_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if hexa20_eids:
                self.write_element_transient('CHEXA',  20,  hexa20_eids, dt, header, hexa_msg, f,is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

            if penta_eids:
                self.write_element_transient('CPENTA', 6, penta_eids, dt, header, penta_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if penta15_eids:
                self.write_element_transient('CPENTA', 15, penta15_eids, dt, header, penta_msg, f,is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
        return pageNum - 1

    def write_element_transient(self, element_name, nnodes, eids, dt, header, msg, f, is_mag_phase):
        dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
        header[1] = dtLine
        msg = header + msg

        f.write('\n'.join(msg))
        for eid in eids:
            node_ids = self.oxx[dt][eid].keys()
            node_ids.remove('CENTER')
            cid = 10
            f.write('0 %12i %11sGRID CS %2i GP\n' % (eid, 0, nnodes))
            for inode in ['CENTER'] + sorted(node_ids):
                # cid
                oxx = self.oxx[dt][eid][inode]
                oyy = self.oyy[dt][eid][inode]
                ozz = self.ozz[dt][eid][inode]
                txy = self.txy[dt][eid][inode]
                tyz = self.tyz[dt][eid][inode]
                txz = self.txz[dt][eid][inode]
                ([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                  oxxi, oyyi, ozzi, txyi, tyzi, txzi,], isAllZeros) = writeImagFloats13E([oxx, oyy, ozz,
                                                                                          txy, tyz, txz], is_mag_phase)

                f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (inode, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('',    oxxi, oyyi, ozzi, txyi, tyzi, txzi))


class ComplexSolidStrainObject(StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)

        self.eType = {}
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.cid = {}  # gridGauss
        self.exx = {}
        self.eyy = {}
        self.ezz = {}
        self.exy = {}
        self.eyz = {}
        self.exz = {}

        #self.dt = dt
        if is_sort1:
            if dt is not None:
                pass
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            raise NotImplementedError('SORT2')
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = []
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.exx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid, exx, eyy, ezz, exy, eyz, exz\n  ')
        msg.append(', '.join(set(self.eType.values()))+'\n  ')
        msg = msg + self.get_data_code()
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            if not hasattr(self, 'data'):
                self.data = []
            self.data += data
        else:
            dt = transient[1]
            if not hasattr(self, 'data'):
                self.data = {}
            #print(self.data)
            if dt not in self.data:
                self.data[dt] = []
            for line in data:
                self.data[dt] += data

    def processF06Data(self):
        raise NotImplementedError()

    def delete_transient(self, dt):
        del self.exx[dt]
        del self.eyy[dt]
        del self.ezz[dt]
        del self.exy[dt]
        del self.eyz[dt]
        del self.exz[dt]

    def get_transients(self):
        k = self.exx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.ezz[dt] = {}
        self.exy[dt] = {}
        self.eyz[dt] = {}
        self.exz[dt] = {}

    def add_eid_sort1(self, element_num, eType, dt, eid, cid, ctype, nodef):
        assert cid >= 0
        assert eid >= 0

        if dt not in self.exx:
            self.add_new_transient(dt)
        assert eid not in self.exx[dt], self.exx[dt]

        self.eType[eid] = eType
        self.cid[eid] = cid
        self.exx[dt][eid] = {}
        self.eyy[dt][eid] = {}
        self.ezz[dt][eid] = {}
        self.exy[dt][eid] = {}
        self.eyz[dt][eid] = {}
        self.exz[dt][eid] = {}

        #msg = "*eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
        #print msg
        #if nodeID == 0:
            #raise ValueError(msg)

    def add_node_sort1(self, dt, eid, node_id, ex, ey, ez, etxy, etyz, etzx):
        #msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, evm)
        assert eid != 'C'
        if node_id == 0:
            node_id = 'CENTER'
        assert node_id != 0, 'CENTER'
        #print "eid=%s nid=%s exx=%s" %(eid, nodeID, exx)
        self.exx[dt][eid][node_id] = ex
        self.eyy[dt][eid][node_id] = ey
        self.ezz[dt][eid][node_id] = ez

        self.exy[dt][eid][node_id] = etxy
        self.eyz[dt][eid][node_id] = etyz
        self.exz[dt][eid][node_id] = etzx

    def getHeaders(self):
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        return headers

    def directionalVectors(self, exx, eyy, ezz, exy, eyz, exz):
        A = [[exx, exy, exz],
             [exy, eyy, eyz],
             [exz, eyz, ezz]]
        (Lambda, v) = eig(A)  # we can't use a hermitian matrix
        return v

    def getF06_Header(self, is_mag_phase):
        if is_mag_phase:
            tetra_msg = [
                '                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (MAGNITUDE/PHASE)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        else:
            tetra_msg = [
                '                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (REAL/IMAGINARY)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        penta_msg = tetra_msg
        hexa_msg = tetra_msg

        tetra_eids = []
        hexa_eids = []
        penta_eids = []

        tetra10_eids = []
        hexa20_eids = []
        penta15_eids = []
        for eid, eType in sorted(self.eType.iteritems()):
            if eType in ['CTETRA4']:
                tetra_eids.append(eid)
            elif eType in ['CTETRA10']:
                tetra10_eids.append(eid)

            elif eType in ['CPENTA']:
                penta_eids.append(eid)
            elif eType in ['CPENTA15']:
                penta_eids.append(eid)

            elif eType in ['CHEXA8']:
                hexa_eids.append(eid)
            elif eType in ['CHEXA20']:
                hexa20_eids.append(eid)

            else:
                raise NotImplementedError('eType=%r' % eType)
        return (tetra_msg, hexa_msg, penta_msg,
                tetra_eids, hexa_eids, penta_eids,
                tetra10_eids, hexa20_eids, penta15_eids)

    def write_f06(self, header, page_stamp, pageNum=1, f=None, is_mag_phase=False):
        (tetra_msg, hexa_msg, penta_msg,
         tetra_eids, hexa_eids, penta_eids,
         tetra10_eids, hexa20_eids, penta15_eids) = self.getF06_Header(is_mag_phase)
        dts = self.exx.keys()
        for dt in sorted(dts):
            if tetra_eids:
                self.write_element_transient('CTETRA', 4, tetra_eids, dt, header, tetra_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if tetra10_eids:
                self.write_element_transient('CTETRA', 10, tetra10_eids, dt, header, tetra_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

            if hexa_eids:
                self.write_element_transient('CHEXA',  8,  hexa_eids, dt, header, hexa_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if hexa20_eids:
                self.write_element_transient('CHEXA',  20,  hexa20_eids, dt, header, hexa_msg, f,is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

            if penta_eids:
                self.write_element_transient('CPENTA', 6, penta_eids, dt, header, penta_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if penta15_eids:
                self.write_element_transient('CPENTA', 15, penta15_eids, dt, header, penta_msg, f,is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

        return pageNum - 1

    def write_element_transient(self, element_name, nnodes, eids, dt, header, msg, f, is_mag_phase):
        dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
        header[1] = dtLine
        msg = header + msg

        f.write('\n'.join(msg))
        for eid in eids:
            node_ids = self.exx[dt][eid].keys()
            node_ids.remove('CENTER')
            cid = 10
            f.write('0 %12i %11sGRID CS %2i GP\n' % (eid, 0, nnodes))
            for inode in ['CENTER'] + sorted(node_ids):
                # cid
                oxx = self.exx[dt][eid][inode]
                oyy = self.eyy[dt][eid][inode]
                ozz = self.ezz[dt][eid][inode]
                txy = self.exy[dt][eid][inode]
                tyz = self.eyz[dt][eid][inode]
                txz = self.exz[dt][eid][inode]
                ([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                  oxxi, oyyi, ozzi, txyi, tyzi, txzi,], isAllZeros) = writeImagFloats13E([oxx, oyy, ozz,
                                                                                          txy, tyz, txz], is_mag_phase)

                f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (inode, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('',    oxxi, oyyi, ozzi, txyi, tyzi, txzi))
