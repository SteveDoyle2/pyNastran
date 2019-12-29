from itertools import cycle
from typing import List
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, OES_Object
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header


class HyperelasticQuadArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        if is_sort1:
            pass
            #if dt is not None:
                #self._add = self._add_sort1
                #self._add_new_eid = self._add_new_eid_sort1
                #self._add_new_node = self._add_new_node_sort1
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    #def get_headers(self):
        #raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def get_headers(self) -> List[str]:
        return ['oxx', 'oyy', 'txy', 'angle', 'majorp', 'minorp']

    #def is_bilinear(self):
        #if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
            #return False
        #elif self.element_type in [144, 64, 82, 70, 75]:  # CQUAD4
            #return True
        #else:
            #raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))

    def build(self):
        """sizes the vectorized attributes of the HyperelasticQuadArray"""
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []

        #if self.element_type in [33, 74]:
            #nnodes_per_element = 1
        #elif self.element_type == 144:
            #nnodes_per_element = 5
        #elif self.element_type == 64:  # CQUAD8
            #nnodes_per_element = 5
        #elif self.element_type == 82:  # CQUADR
            #nnodes_per_element = 5
        #elif self.element_type == 70:  # CTRIAR
            #nnodes_per_element = 4
        #elif self.element_type == 75:  # CTRIA6
            #nnodes_per_element = 4
        #else:
        if self.element_type == 139: # QUAD4FD
            nnodes_per_element = 4
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

        #print('nnodes_per_element[%s, %s] = %s' % (self.isubcase, self.element_type, nnodes_per_element))
        #self.nnodes = nnodes_per_element
        #self.nelements //= nnodes_per_element
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.element_node = np.zeros((self.ntotal, 2), dtype='int32')

        #self.Type[eid] = Type
        #self.oxx[dt] = {eid: [oxx]}
        #self.oyy[dt] = {eid: [oyy]}
        #self.txy[dt] = {eid: [txy]}
        #self.angle[dt] = {eid: [angle]}
        #self.majorP[dt] = {eid: [majorP]}
        #self.minorP[dt] = {eid: [minorP]}


        #[oxx, oyy, txy, angle, majorp, minorp]
        self.data = np.zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    #def build_dataframe(self):
        #"""creates a pandas dataframe"""
        #import pandas as pd
        #headers = self.get_headers()

        #nelements = self.element_node.shape[0] // 2
        #if self.is_fiber_distance:
            #fiber_distance = ['Top', 'Bottom'] * nelements
        #else:
            #fiber_distance = ['Mean', 'Curvature'] * nelements
        #fd = np.array(fiber_distance, dtype='unicode')
        #element_node = [self.element_node[:, 0], self.element_node[:, 1], fd]

        #if self.nonlinear_factor not in (None, np.nan):
            #column_names, column_values = self._build_dataframe_transient_header()
            #self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
            #self.data_frame.index.names = ['ElementID', 'NodeID', 'Location', 'Item']
        #else:
            ## option B - nice!
            #df1 = pd.DataFrame(element_node).T
            #df1.columns = ['ElementID', 'NodeID', 'Location']
            #df2 = pd.DataFrame(self.data[0])
            #df2.columns = headers
            #self.data_frame = df1.join(df2)
        #self.data_frame = self.data_frame.reset_index().replace({'NodeID': {0:'CEN'}}).set_index(['ElementID', 'NodeID', 'Location'])
        #print(self.data_frame)

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fiber_dist1, oxx1, oyy1, txy1, angle1, majorP1, minorP1, ovm1) = t1
                    (fiber_dist2, oxx2, oyy2, txy2, angle2, majorP2, minorP2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            fiber_dist1, oxx1, oyy1, txy1, angle1, majorP1, minorP1, ovm1,
                            fiber_dist2, oxx2, oyy2, txy2, angle2, majorP2, minorP2, ovm2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    #def _add_new_eid(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #self._add_new_eid_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)


    def _add_new_eid_sort1(self, dt, eid, eype, node_id, oxx, oyy, txy, angle, majorP, minorP):
        assert isinstance(eid, integer_types), eid
        assert isinstance(node_id, integer_types), node_id
        self._times[self.itime] = dt
        #assert self.itotal == 0, oxx
        self.element_node[self.itotal, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [oxx, oyy, txy, angle, majorP, minorP]
        self.itotal += 1
        self.ielement += 1

    #def _add_new_node(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert isinstance(node_id, integer_types), node_id
        #self._add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    #def _add(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert isinstance(node_id, integer_types), node_id
        #self._add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    #def _add_new_node_sort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #self._add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm)

    def _add_sort1(self, dt, eid, etype, node_id, oxx, oyy, txy, angle, majorP, minorP):
        assert eid is not None, eid
        assert isinstance(node_id, integer_types), node_id
        self.element_node[self.itotal, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [oxx, oyy, txy, angle, majorP, minorP]
        self.itotal += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = 4
        ntotal = self.ntotal
        nlayers = 2
        nelements = self.ntotal // nnodes

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msgi = '  type=%s ntimes=%i nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, nlayers, ntotal)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, nelements, nnodes, nlayers, ntotal)
            ntimes_word = '1'
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 str(', '.join(headers))))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  data.shape=%s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg.append('  s_code: %s\n' % self.s_code)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = np.searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = np.ravel([np.searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #msg, nnodes, cen = _get_plate_msg(self)


        msg = ['           S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUAD4FD )\n',
               '  ELEMENT     GRID/    POINT       ---------CAUCHY STRESSES--------             PRINCIPAL STRESSES (ZERO SHEAR)\n',
               '     ID       GAUSS      ID      NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR\n', ]
            #0       1     GAUS         1   7.318995E+00   6.367099E-01  -6.551054E+00   -31.4888    1.133173E+01   -3.376026E+00
            #                           2   1.097933E+01   4.149028E+00   6.278160E+00    30.7275    1.471111E+01    4.172537E-01

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[oxx, oyy, txy, angle, majorp, minorp]
            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            txy = self.data[itime, :, 2]
            angle = self.data[itime, :, 3]
            majorP = self.data[itime, :, 4]
            minorP = self.data[itime, :, 5]

            for (i, eid, nid, oxxi, oyyi, txyi, anglei, major, minor) in zip(
                 cycle([1, 2, 3, 4]), eids, nids, oxx, oyy, txy, angle, majorP, minorP):
                [oxxi, oyyi, txyi, major, minor] = write_floats_13e(
                    [oxxi, oyyi, txyi, major, minor])

                if i == 1:
                    gauss = 'GAUS'  # TODO: update
                    f06_file.write(
                        '0%8i %8s  %8i  %13s  %13s  %13s  %13s  %13s  %13s\n' % (
                            eid, gauss, 1, oxxi, oyyi, txyi, anglei, major, minor))
                else:
                    f06_file.write(
                        ' %8s %8s  %8i  %13s  %13s  %13s  %13s  %13s  %13s\n' % (
                            '', '', i, oxxi, oyyi, txyi, anglei, major, minor))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    #def get_nnodes_bilinear(self):
        #is_bilinear = False
        #if self.element_type == 74:
            #nnodes = 3
        #elif self.element_type == 33:
            #nnodes = 4
        #elif self.element_type == 144:
            #nnodes = 4
            #is_bilinear = True
        #elif self.element_type == 82:  # CQUADR
            #nnodes = 4
            #is_bilinear = True
        #elif self.element_type == 64:  # CQUAD8
            #nnodes = 4
            #is_bilinear = True
        #elif self.element_type == 75:  # CTRIA6
            #nnodes = 3
            #is_bilinear = True
        #elif self.element_type == 70:  # CTRIAR
            #nnodes = 3
            #is_bilinear = True
        #else:
            #raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))
        #return nnodes, is_bilinear
