# pylint: disable=C0301,C0103,R0913,R0914,R0904,C0111,R0201,R0902
import warnings
from itertools import count
from struct import pack

import numpy as np
from numpy import zeros

#from pyNastran.utils.numpy_utils import float_types
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.result_objects.utils_pandas import build_dataframe_transient_header, build_pandas_transient_element_node
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.op2.op2_interface.write_utils import to_column_bytes


class RealSolidCompositeArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            ##sort1
            #self.add_node = self.add_node_sort1
            #self.add_eid = self.add_eid_sort1
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    # def update_data_components(self):
    #     ntimes, nelements_nnodes = self.data.shape[:2]
    #     # vm
    #     oxx = self.data[:, :, 0].reshape(ntimes * nelements_nnodes)
    #     oyy = self.data[:, :, 1].reshape(ntimes * nelements_nnodes)
    #     ozz = self.data[:, :, 2].reshape(ntimes * nelements_nnodes)
    #     txy = self.data[:, :, 3].reshape(ntimes * nelements_nnodes)
    #     tyz = self.data[:, :, 4].reshape(ntimes * nelements_nnodes)
    #     txz = self.data[:, :, 5].reshape(ntimes * nelements_nnodes)
    #
    #     #I1 = oxx + oyy + ozz
    #     #txyz = txy**2 + tyz**2 + txz ** 2
    #     #I2 = oxx * oyy + oyy * ozz + ozz * oxx - txyz
    #     #I3 = oxx * oyy * ozz + 2 * txy * tyz * txz + oxx * tyz**2 - oyy * txz**2 - ozz * txy
    #
    #     # (n_subarrays, nrows, ncols)
    #     o1, o2, o3 = principal_components_3d(
    #         ntimes, nelements_nnodes,
    #         oxx, oyy, ozz, txy, tyz, txz,
    #         self.is_stress)
    #     ovm_sheari = ovm_shear_3d(oxx, oyy, ozz, txy, tyz, txz, o1, o3,
    #                               self.is_von_mises, self.is_stress)
    #     ovm_sheari2 = ovm_sheari.reshape(ntimes, nelements_nnodes)
    #
    #     self.data[:, :, 6] = o1.reshape(ntimes, nelements_nnodes)
    #     self.data[:, :, 7] = o2.reshape(ntimes, nelements_nnodes)
    #     self.data[:, :, 8] = o3.reshape(ntimes, nelements_nnodes)
    #     self.data[:, :, 9] = ovm_sheari2
    #
    #     A = [[doxx, dtxy, dtxz],
    #          [dtxy, doyy, dtyz],
    #          [dtxz, dtyz, dozz]]
    #     (_lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

    def build(self):
        """sizes the vectorized attributes of the RealSolidArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
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

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)

        if self.is_sort1:
            ntimes = self.ntimes
            ntotal = self.ntotal
            nelements = self.nelements
        else:
            #print(f'ntimes={self.ntimes} nelements={self.nelements} ntotal={self.ntotal}')
            ntimes = self.nelements
            ntotal = self.ntotal
            nelements = self.ntimes
            #print(f'ntimes={ntimes} nelements={nelements} ntotal={ntotal}')
        #self.ntimes = ntimes
        #self.ntotal = ntotal
        #self.nelements = nelements

        _times = zeros(ntimes, dtype=self.analysis_fmt)

        # TODO: could be more efficient by using nelements for cid
        #element_cid = zeros((nelements, 2), dtype=idtype)
        #if nelements > 5000:
            #raise RuntimeError(nelements)

        #if self.element_name == 'CTETRA':
            #nnodes = 4
        #elif self.element_name == 'CPENTA':
            #nnodes = 6
        #elif self.element_name == 'CHEXA':
            #nnodes = 8
        #self.element_node = zeros((self.ntotal, nnodes, 2), 'int32')

        #eid_device, layer,
        #[o1, o2, t12, t1z, t2z, angle, major, minor, ovm]
        if self.num_wide == 11:
            nodes_per_element = 1
        elif self.num_wide == 43:
            nodes_per_element = 5
        else:
            raise NotImplementedError(self.code_information())

        ntotal *= nodes_per_element

        element_layer_node = np.zeros((ntotal, 3), dtype=idtype)
        data = np.zeros((ntimes, ntotal, 7), fdtype)
        #print('RealSolidCompositeArray: data.shape=%s' % str(data.shape))
        #self.nnodes = element_layer.shape[0] // self.nelements
        #self.data = zeros((self.ntimes, self.nelements, nnodes+1, 10), 'float32')

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element_layer_node = group.create_dataset('element_layer_node', data=element_layer_node)
            #self.element_cid = group.create_dataset('element_cid', data=element_cid)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element_layer_node = element_layer_node
            #self.element_cid = element_cid
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.headers
        # TODO: cid?
        #element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = build_dataframe_transient_header(self)
            data_frame = build_pandas_transient_element_node(
                self, column_values, column_names,
                headers, self.element_layer_node, self.data)
            #self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
            #self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']
            print(data_frame)
            raise RuntimeError('finish pd.Panel')
        else:
            # Static                     oxx       oyy  ...       t2z  von_mises
            # ElementID Layer Grid                      ...
            # 1         1     0     1.060943  0.179359  ...  0.019784   0.885374
            #           2     0     1.089592  0.183899  ...  0.083607   0.921374
            # 2         1     0     1.107400  0.044122  ... -0.011024   1.071296
            #           2     0     1.043136  0.033165  ... -0.038763   1.012717
            index = pd.MultiIndex.from_arrays(self.element_layer_node.T, names=['ElementID', 'Layer', 'Grid'])
            data_frame = pd.DataFrame(self.data[0, :, :], columns=headers, index=index)
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def add_sort1(self, dt, eid, layer, location, grid, o11, o22, o33, t12, t2z, t1z, ovm):
        assert self.sort_method == 1, self
        # See the CHEXA, CPENTA, or CTETRA entry for the definition of the element coordinate systems.
        # The material coordinate system (CORDM) may be the basic system (0 or blank), any defined system
        # (Integer > 0), or the standard internal coordinate system of the element designated as:
        # -1: element coordinate system (-1)
        # -2: element system based on eigenvalue techniques to insure non bias in the element formulation(-2).
        #     C:\MSC.Software\msc_nastran_runs\ecs-2-rg.op2
        #assert cid >= -2, cid
        assert eid >= 0, eid
        assert layer >= 1, layer

        #print(f'dt={dt} eid={eid}')
        self._times[self.itime] = dt
        self.element_layer_node[self.itotal, :] = [eid, layer, grid]  # 0 is center

        self.data[self.itime, self.itotal, :] = [o11, o22, o33, t12, t2z, t1z, ovm]

        #print('element_cid[%i, :] = [%s, %s]' % (self.ielement, eid, cid))
        #if self.ielement == self.nelements:
            #self.ielement = 0
        #self.element_cid[self.ielement, :] = [eid, cid]
        self.itotal += 1
        self.ielement += 1

    #def add_node_sort1(self, dt, eid, unused_inode, node_id,
                       #oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                       #unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        ## skipping aCos, bCos, cCos, pressure
        #omax_mid_min = [o1, o2, o3]
        #omin = min(omax_mid_min)
        #omax_mid_min.remove(omin)

        #omax = max(omax_mid_min)
        #omax_mid_min.remove(omax)

        #omid = omax_mid_min[0]
        #self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]
        ##print('data[%s, %s, :] = %s' % (self.itime, self.itotal, str(self.data[self.itime, self.itotal, :])))

        ##self.data[self.itime, self.ielement-1, self.inode, :] = [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]

        ##print('eid=%i node_id=%i exx=%s' % (eid, node_id, str(oxx)))
        #self.element_node[self.itotal, :] = [eid, node_id]
        ##self.element_node[self.ielement-1, self.inode-1, :] = [eid, node_id]
        #self.itotal += 1

    #def add_eid_sort2(self, unused_etype, cid, dt, eid, unused_node_id,
                      #oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                      #unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        ##itime = self.ielement
        ##ielement = self.itotal
        ##itotal = self.itime
        ##print(self.ntimes, self.nelements, self.ntotal, self.nnodes)
        #itime = self.itotal // self.nnodes
        #ielement = self.itime
        #itotal = self.itotal
        #assert cid >= -2, cid
        #assert eid >= 0, eid

        ##try:
        #self._times[itime] = dt
            ##print(f'dt={dt} eid={eid} ielement={ielement} -> itime={itime} itotal={itotal}')
        ##except IndexError:
            ##print(f'*dt={dt} eid={eid} ielement={ielement} -> itime={itime} itotal={itotal}')
            ##self.itime += 1
            ##self.ielement += 1
            ##return
        #self.element_node[itotal, :] = [eid, 0]  # 0 is center

        #omax_mid_min = [o1, o2, o3]
        #omin = min(omax_mid_min)
        #omax_mid_min.remove(omin)

        #omax = max(omax_mid_min)
        #omax_mid_min.remove(omax)

        #omid = omax_mid_min[0]
        #self.data[itime, itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]

        ##print('element_cid[%i, :] = [%s, %s]' % (self.ielement, eid, cid))
        ##if self.ielement == self.nelements:
            ##self.ielement = 0
        #self.element_cid[ielement, :] = [eid, cid]
        ##self.itime += 1
        #self.itotal += 1
        #self.ielement += 1
        ##print('self._times', self._times)

    #def add_node_sort2(self, dt, eid, unused_inode, node_id,
                       #oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                       #unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        ##ielement = self.ielement
        ##itotal = self.itotal
        ##itime = self.itime
        ##itime=0 ielement=1 itotal=1
        ##itime=0 ielement=1 itotal=2
        ##itime=0 ielement=1 itotal=3
        ##itime=0 ielement=1 itotal=4

        ##ielement = self.ielement
        ##itime = (self.itime - 1) % self.nelements
        ##itime = self.itime - 1
        #nnodes = self.nnodes
        #itime = self.itotal // nnodes
        #itotal = self.itotal
        ##ielement = self.ielement - 1
        ##ielement = self.itime
        ##inode = self.itotal % nnodes
        ##itotal2 = (self.ielement - 1) * nnodes + inode
        ##print(f'  itime={itime} itotal={itotal}; nid={node_id}; '
              ##f'ielement={ielement} inode={inode} -> itotal2={itotal2}')

        ## skipping aCos, bCos, cCos, pressure
        #omax_mid_min = [o1, o2, o3]
        #omin = min(omax_mid_min)
        #omax_mid_min.remove(omin)

        #omax = max(omax_mid_min)
        #omax_mid_min.remove(omax)

        #omid = omax_mid_min[0]
        #self.data[itime, itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]
        ##print('data[%s, %s, :] = %s' % (self.itime, self.itotal, str(self.data[self.itime, self.itotal, :])))

        ##print('eid=%i node_id=%i exx=%s' % (eid, node_id, str(oxx)))
        #self.element_node[itotal, :] = [eid, node_id]
        ##self.element_node[ielement-1, inode-1, :] = [eid, node_id]
        #self.itotal += 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ieid, eid_layer_nid in enumerate(self.element_layer_node):
                    (eid, layer, nid) = eid_layer_nid
                    t1 = self.data[itime, ieid, :]
                    t2 = table.data[itime, ieid, :]
                    (oxx1, oyy1, ozz1, txy1, tyz1, txz1, o11, o21, o31, ovm1) = t1
                    (oxx2, oyy2, ozz2, txy2, tyz2, txz2, o12, o22, o32, ovm2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '%s      (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid, nid,
                                oxx1, oyy1, ozz1, txy1, tyz1, txz1, o11, o21, o31, ovm1,
                                ' ' * (len(str(eid)) + len(str(nid)) + 2),
                                oxx2, oyy2, ozz2, txy2, tyz2, txz2, o12, o22, o32, ovm2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    @property
    def nnodes_per_element(self) -> int:
        return self.nnodes_per_element_no_centroid + 1

    @property
    def nnodes_per_element_no_centroid(self) -> int:
        if self.element_type == 269: # CHEXA
            nnodes = 8
        elif self.element_type == 270: # CPENTA
            nnodes = 6
        else:
            raise NotImplementedError(f'element_name={self.element_name} self.element_type={self.element_type}')
        return nnodes

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        #nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        #try:
            #nnodes_per_element = self.element_layer_node.shape[0] // nelements
        #except ZeroDivisionError:
            #nnodes_per_element = '???'
        eids = self.element_layer_node[:, 0]
        layers = self.element_layer_node[:, 1]

        nelements = len(np.unique(eids))
        #nids = self.element_layer_node[:10, 2]
        eid0 = eids[0]
        ieid = np.where(eid0 == eids[:10])[0]
        ilayer0 = layers[ieid[-1]]
        nnodes = len(ieid) // ilayer0

        nlayers = len(eids)
        nlayers_max = layers.max()
        #nnodes = self.element_layer_grid.shape[0]


        msg = []

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%d nelements=%d nlayers=%d nlayers_max=%d\n  nnodes_per_element=%d (including centroid)\n'
                       % (self.__class__.__name__, ntimes, nelements, nlayers, nlayers_max, nnodes))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%d nlayers=%d nlayers_max=%d\n  nodes_per_element=%d (including centroid)\n'
                       % (self.__class__.__name__, nelements, nlayers, nlayers_max, nnodes))
            ntimes_word = '1'
        msg.append('  eType, cid\n')
        headers = self.headers
        n = len(headers)
        msg.append('  data: [%s, nlayers, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  element_layer_node.shape = {self.element_layer_node.shape}\n')
        #msg.append(f'  element_cid.shape = {self.element_cid.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        #print(''.join(msg))
        return msg


    #def get_element_index(self, eids):
        ## elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        ##ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element_node[:, 0])
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    def write_f06(self, f06_file, header=None, page_stamp: str='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        nnodes, msg_temp = _get_f06_header_nnodes(self, is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_layer_node[:, 0]
        layers = self.element_layer_node[:, 1]
        nodes = self.element_layer_node[:, 2]

        #eids3 = self.element_cid[:, 0]
        #cids3 = self.element_cid[:, 1]

        fdtype = self.data.dtype
        #o11, o22, o33, t12, t23, t13, ovm
        o11 = self.data[:, :, 0]
        o22 = self.data[:, :, 1]
        o33 = self.data[:, :, 2]
        t12 = self.data[:, :, 3]
        t23 = self.data[:, :, 4]
        t13 = self.data[:, :, 5]
        ovm = self.data[:, :, 6]

        #nnodes_total = self.data.shape[1]
        #print(self.element_layer_node[:, 0].tolist())
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            o11 = self.data[itime, :, 0]
            o22 = self.data[itime, :, 1]
            o33 = self.data[itime, :, 2]
            t12 = self.data[itime, :, 3]
            t23 = self.data[itime, :, 4]
            t13 = self.data[itime, :, 5]
            ovm = self.data[itime, :, 6]

            #cnnodes = nnodes + 1
            for i, deid, layer, node_id, do11, do22, do33, dt12, dt23, dt13, dovm in zip(
                    count(), eids, layers, nodes, o11, o22, o33, t12, t23, t13, ovm):

                # o1-max
                # o2-mid
                # o3-min
                # 1.060943E+00
                # -2.357560E-03
                [o11i, o22i, o33i, t12i, t23i, t13i, ovmi] = write_floats_13e(
                    [do11, do22, do33, dt12, dt23, dt13, dovm])

                if node_id == 0:
                    f06_file.write('0  %8d  %6d   MID   CENTER %s %s %s %s %s %s %s\n' % (deid, layer, o11i, o22i, o33i, t12i, t23i, t13i, ovmi))
                else:
                    f06_file.write('                          %8d %s %s %s %s %s %s %s\n' % (node_id, o11i, o22i, o33i, t12i, t23i, t13i, ovmi))
                i += 1
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def __write_op2(self, op2_file, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        calculate_directional_vectors = True
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            #print('***************', itable)
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)
        nnodes_expected = self.nnodes

        eids2 = self.element_layer_node[:, 0]
        layer = self.element_layer_node[:, 1]
        nodes = self.element_layer_node[:, 2]
        nelements_nodes = len(nodes)

        #eids3 = self.element_cid[:, 0]
        #cids3 = self.element_cid[:, 1]
        #element_device = eids3 * 10 + self.device_code

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nelements = len(np.unique(eids2))

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)
        nnodes_centroid = self.nnodes_per_element
        nnodes_no_centroid = self.nnodes_per_element_no_centroid
        ntotali = 4 + 21 * nnodes_no_centroid
        ntotali = self.num_wide
        ntotal = ntotali * nelements


        #print('shape = %s' % str(self.data.shape))
        assert nnodes > 1, nnodes
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write(f'  ntimes = {self.ntimes}\n')
        ntimes = self.ntimes

        #print('ntotal=%s' % (ntotal))
        if not self.is_sort1:
            raise NotImplementedError('SORT2')
        #op2_format = endian + b'2i6f'

        idtype = self.element_cid.dtype
        fdtype = self.data.dtype
        if self.size == fdtype.itemsize:
            grid_bytes = b'GRID'
        else:
            # warnings.warn(f'downcasting {self.class_name}...')
            idtype = np.int32(1)
            fdtype = np.float32(1.0)
            grid_bytes = b'GRID'

        cen_array = np.full(nelements, grid_bytes, dtype='|S4')
        nnodes_no_centroid_array = np.full(nelements, nnodes_no_centroid, dtype=idtype)

        element_wise_data = to_column_bytes([
            element_device, # ints
            cids3, # ints
            cen_array, # bytes
            nnodes_no_centroid_array, # ints
        ], fdtype, debug=False)

        oxx = self.data[:, :, 0]
        oyy = self.data[:, :, 1]
        ozz = self.data[:, :, 2]
        txy = self.data[:, :, 3]
        tyz = self.data[:, :, 4]
        txz = self.data[:, :, 5]
        o1 = self.data[:, :, 6]
        o2 = self.data[:, :, 7]
        o3 = self.data[:, :, 8]
        ovm = self.data[:, :, 9]
        p = (o1 + o2 + o3) / -3.

        # speed up transient cases, but slightly slows down static cases
        data_out = np.empty((nelements, 4+21*nnodes_centroid), dtype=fdtype)

        # setting:
        #  - CTETRA: [element_device, cid, 'CEN/', 4]
        #  - CPYRAM: [element_device, cid, 'CEN/', 5]
        #  - CPENTA: [element_device, cid, 'CEN/', 6]
        #  - CHEXA:  [element_device, cid, 'CEN/', 8]
        data_out[:, :4] = element_wise_data

        # we could tack the nodes on, so we don't have to keep stacking it
        # but we run into issues with datai
        #
        # total=nelements_nodes
        #nodes_view = nodes.view(fdtype).reshape(nelements, nnodes_centroid)
        #inode = np.arange(nnodes_centroid)
        #data_out[:, 4+inode*21] = nodes_view[:, inode]

        # v is the (3, 3) eigenvector for every time and every element
        if calculate_directional_vectors:
            v = calculate_principal_eigenvectors4(
                ntimes, nnodes,
                oxx, oyy, ozz, txy, txz, tyz,
                fdtype)[1]
        else:
            v = np.zeros((ntimes, nnodes, 3, 3), dtype=fdtype)

        op2_ascii.write(f'nelements={nelements:d}\n')
        for itime in range(self.ntimes):
            vi = v[itime, :, :, :]
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            col_inputs = [
                nodes,
                oxx[itime, :], txy[itime, :], o1[itime, :], vi[:, 0, 1], vi[:, 0, 2], vi[:, 0, 0], p[itime, :], ovm[itime, :],
                oyy[itime, :], tyz[itime, :], o2[itime, :], vi[:, 1, 1], vi[:, 1, 2], vi[:, 1, 0],
                ozz[itime, :], txz[itime, :], o3[itime, :], vi[:, 2, 1], vi[:, 2, 2], vi[:, 2, 0],
            ]

            # stack each output by columns and fix any dtypes
            datai = to_column_bytes(col_inputs, fdtype)
            #datai2 = datai.reshape(nelements, 21*nnodes_centroid)
            #data_out = np.hstack([element_wise_data, datai2])
            #data_out[:, 4:] = datai2

            # switch datai to element format and put it in the output buffer
            data_out[:, 4:] = datai.reshape(nelements, 21*nnodes_centroid)
            op2_file.write(data_out)

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealSolidCompositeStressArray(RealSolidCompositeArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidCompositeArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    @property
    def headers(self) -> list[str]:
        if self.is_von_mises:
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['oxx', 'oyy', 'ozz', 't12', 't1z', 't2z', von_mises]
        return headers


class RealSolidCompositeStrainArray(RealSolidCompositeArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidCompositeArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    @property
    def headers(self) -> list[str]:
        if self.is_von_mises:
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz', von_mises]
        return headers

def _get_solid_msgs(self):
    """
    '                      S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( H E X A )'
    ' '
    '0   ELEMENT     PLY FIBER    EDGE                                                                                                    '
    '      ID         ID   LOC  GRID ID   NORMAL-1      NORMAL-2      NORMAL-3      SHEAR-12      SHEAR-23      SHEAR-13      VON MISES   '
    '0         1       1   MID   CENTER  1.060943E+00  1.793585E-01  1.853168E-01  5.978239E-02 -2.357560E-03  1.978396E-02  8.853740E-01'
    """
    if self.is_von_mises:
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    tetra_msg = None
    penta_msg = None
    pyram_msg = None
    if self.is_stress:
        base_msg = [
            ' \n'
            '0   ELEMENT     PLY FIBER    EDGE                                                                                                    \n'
            '      ID         ID   LOC  GRID ID   NORMAL-1      NORMAL-2      NORMAL-3      SHEAR-12      SHEAR-23      SHEAR-13      %s   \n' % von_mises
        ]
        #tetra_msg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        #penta_msg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                      S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( H E X A )\n', ]
        #pyram_msg = ['                      S T R E S S E S   I N   P Y R A M I D   S O L I D   E L E M E N T S   ( P Y R A M )\n', ]
    else:
        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT  STRAINS---------       DIR.  COSINES       MEAN                   \n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises]
        #tetra_msg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        #penta_msg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                        S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( H E X A )\n', ]
        #pyram_msg = ['                        S T R A I N S   I N   P Y R A M I D   S O L I D   E L E M E N T S   ( P Y R A M )\n', ]

    #tetra_msg += base_msg
    #penta_msg += base_msg
    hexa_msg += base_msg
    return tetra_msg, penta_msg, hexa_msg, pyram_msg

def _get_f06_header_nnodes(self, is_mag_phase=True):
    tetra_msg, penta_msg, hexa_msg, pyram_msg = _get_solid_msgs(self)
    if self.element_type == 269: # CHEXA
        msg = hexa_msg
        nnodes = 8
    elif self.element_type == 270: # CPENTA
        msg = penta_msg
        nnodes = 6
    else:  # pragma: no cover
        msg = f'element_name={self.element_name} self.element_type={self.element_type}'
        raise NotImplementedError(msg)
    return nnodes, msg
